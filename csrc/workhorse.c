#include <Python.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

#ifdef _WIN32 // 32/64 bit windows
#pragma message ( "*** compiling for windows ***" )
#include <windows.h>
#endif

#define SCANBUFSIZE (1024*1024)
#define MIN(A,B) (A<B?A:B)
#define MAX(A,B) (A>B?A:B)

#if 0
#define DBG(fmt, ...) fprintf(stderr, fmt "\n", __VA_ARGS__)
#else
#define DBG(...)
#endif

int running=0, stop=0;
// see engine_config()
int maxerrors=0;
int minoverlap=20;
int minreadlength=10;
int nthreads=1;
char Amin='!', Azero='!';
int domore=0;
// synchronizing
pthread_mutex_t ll_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t rl_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t args_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t records_parsed_mutex = PTHREAD_MUTEX_INITIALIZER;
// python objects for interfacing etc
PyObject *engine_mod, *hittuple;
// python pthreads exception handling
PyObject *exception, *fastq_exception;
#define ERRSTR_LENGTH 1024
char errstr[ERRSTR_LENGTH];

// stats
long total, parsed, sigints;
int nseqs;
long *seqbasehits, *seqhits, records_parsed;
long nrecords;
long nN, nG, nA, nC, nT, nX;
#define MAX_READLENGTH 1024
long rls_longest;
long rls_buf[MAX_READLENGTH];

void add_rl(long rl)
{
    pthread_mutex_lock(&rl_mutex);
    if (rl>=0 && rl<MAX_READLENGTH)
        rls_buf[rl]++;
    if (rl > rls_longest)
        rls_longest = rl;
    pthread_mutex_unlock(&rl_mutex);
}

int Amin_steps = 5;
long **all_rls_buf;
void analyse_record(char *rstart, long blen)
{
    int i;
    char *ptr,*qtr;
    nrecords++;
    for(ptr=rstart, i=0; i<4 && ptr-rstart<blen; ptr++) {
        if (*ptr == '\n')
            i++;
        if (i == 0) switch(*ptr) {
            case 'N': nN += 1; break;
            case 'A': nA += 1; break;
            case 'G': nG += 1; break;
            case 'T': nT += 1; break;
            case 'C': nC += 1; break;
            default: nX += 1; break;
        }
    }
    if (i<4) return;
    for(ptr=rstart; i<3; ptr++)
        if (*ptr == '\n')
            i++;
    for(i=0; i<2*Amin_steps; i++)
        for(++ptr, qtr=ptr; *ptr!='\n'; ptr++)
            if (*ptr<(Amin+(i<Amin_steps?-i-1:i-Amin_steps+1))) {
                if (qtr)
                    all_rls_buf[i][qtr-ptr>=MAX_READLENGTH ? MAX_READLENGTH-1 : ptr-qtr]++;
                qtr = NULL;
            } else if (!qtr) qtr=ptr;
}

void init_stats(char **seqlist)
{
    int i;

    //FIXME : allocated structures are never freed

    memset((void *) rls_buf, 0, sizeof(rls_buf));
    rls_longest = -1;

    nrecords = nA = nG = nC = nT = nN = nX = 0;

    if (all_rls_buf) {
        for(i=0; i<2*Amin_steps; i++)
            free(all_rls_buf[i]);
        free(all_rls_buf);
    }
    all_rls_buf = (long **) malloc(sizeof(long *) * Amin_steps*2);
    for(i=0; i<2*Amin_steps; i++) {
        all_rls_buf[i] = (long *) malloc(sizeof(long) * MAX_READLENGTH);
        memset((void *) all_rls_buf[i], 0, sizeof(long) * MAX_READLENGTH);
    }

    for(nseqs=0; seqlist[nseqs]; nseqs++);
    seqbasehits = (long *) malloc(sizeof(long) * nseqs);
    memset((void *) seqbasehits, 0, sizeof(long) * nseqs);
    seqhits = (long *) malloc(sizeof(long) * nseqs);
    memset((void *) seqhits, 0, sizeof(long) * nseqs);
    records_parsed = 0;
}

static PyObject *
engine_stats(PyObject *self, PyObject *args)
{
    int i;
    PyObject *rls, *sbhs, *shs;

    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    rls = PyTuple_New(rls_longest+1);
    for(i=0; i<=rls_longest; i++) {
        PyTuple_SetItem(rls, i, PyInt_FromLong(rls_buf[i]));
        //DBG("rl=%d %ldx", i, rls_buf[i]);
    }

    sbhs = PyTuple_New(nseqs);
    for(i=0; i<nseqs; i++)
        PyTuple_SetItem(sbhs, i, PyInt_FromLong(seqbasehits[i]));

    shs = PyTuple_New(nseqs);
    for(i=0; i<nseqs; i++)
        PyTuple_SetItem(shs, i, PyInt_FromLong(seqhits[i]));

    // DBG("engine_stats : parsed=%li total=%li", parsed, total);

    return Py_BuildValue("{sOsfsOsOsOsOsOsO}",
            "readlengths", rls,
            "progress", ((float) (parsed>total?total:parsed))/(total ? total : 1),
            "nseqbasehits", sbhs,
            "nseqhits", shs,
            "parsed", PyInt_FromLong(parsed),
            "total", PyInt_FromLong(total),
            "sigints", PyInt_FromLong(sigints),
            "records_parsed", PyInt_FromLong(records_parsed)
        );
}

void add_records_parsed(long n) {
    pthread_mutex_lock(&records_parsed_mutex);
    DBG("records_parsed : %ld -> %ld", records_parsed, records_parsed + n);
    records_parsed += n;
    pthread_mutex_unlock(&records_parsed_mutex);
}

// hits
typedef struct ll_item_ {
    int seqi;
    long fpos;
    int spos;
    int length;
    int readlength;
    struct ll_item_ *next;
} ll_item;

ll_item *add_ll(ll_item *root, int seqi, long fpos, int spos, int length, int readlength)
{
    pthread_mutex_lock(&ll_mutex);

    //DBG("adding item seqi=%d fpos=%li spos=%i length=%i (thread %li)",
    //        seqi, fpos, spos, length, thread_self());

    while(root->next != NULL)
        root = root->next;
    root->next = (ll_item *) malloc(sizeof(ll_item));

    if (root->next == NULL)
    {
        pthread_mutex_unlock(&ll_mutex);
        return (ll_item *) PyErr_NoMemory();
    }

    root->next->seqi = seqi;
    root->next->fpos = fpos;
    root->next->spos = spos;
    root->next->length = length;
    root->next->readlength = readlength;
    root->next->next = NULL;
    // DBG("added %ld at %ld (thread %ld)", which, pos, thread_self());

    seqbasehits[seqi] += length;
    seqhits[seqi]++;

    pthread_mutex_unlock(&ll_mutex);
    return root->next;
}

void free_ll(ll_item *root)
{
    ll_item *tmp;

    while(root !=NULL)
    {
        tmp = root->next;
        free(root);
        root = tmp;
    }
}


#ifdef _WIN32 // 32/64 bit windows
BOOL CtrlHandler(DWORD fdwCtrlType)
{
    DBG("CtrlHandler called");
    if (fdwCtrlType == CTRL_C_EVENT)
    {
        DBG("CtrlHandler detected CTRL_C_EVENT");
        sigints++;
        return TRUE;
    } 
    return FALSE;
} 
long thread_self() {
    return (long) (pthread_self()).p;
}
long pthread2long(pthread_t x) {
    return (long) (pthread_self()).p;
}
#else // POSIX
void sigint_cb(int sig)
{
    sigints++;
}
long thread_self() {
    return (long) pthread_self();
}
long pthread2long(pthread_t x) {
    return (long) x;
}
#endif

static PyObject *
engine_stop(PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    stop++;

    Py_RETURN_NONE;
}


struct scanargs {
    int assigned;
    char *fname; // NULL indicates end of array
    long start;
    long length;
    char **seqlist;
    int *seqlengths;
    ll_item *root;
    int hits;
};

void dump_record(char *startrecord, char *startread, int rl) {
    int i, j;
    char t, *bases;

    for(i=0,j=0,bases=NULL; startrecord[i] && j<4; i++)
        if (startrecord[i] == '\n') {
            j++;
            if (j == 1)
                bases = startrecord+i+1;
        }
    if (j<4) {
        fprintf(stderr,"*** cannot dump record; j=%d ***\n", j);
        return;
    }

    fprintf(stderr,"<<<<<<<<<<<<<<<<\n");

    t = startrecord[i];
    startrecord[i] = 0;
    fprintf(stderr,"%s", startrecord);
    startrecord[i] = t;

    for(i=0; i<startread-bases; i++)
        fputc(' ', stderr);
    for(i=0; i<rl; i++)
        fputc('*', stderr);
    fputc('\n', stderr);

    fprintf(stderr,">>>>>>>>>>>>>>>>\n");
}

/** 
 * returns the minimal number of bytes that have to be discarded at the beginning
 * of buf to start with the beginning of a fastq record
 * @param buf character buffer containing fastq records
 * @param bsize length of data in buf
 * @return location of beginning of first complete record in buffer -- -1 if no
 *         record could be found
 */

long forward_fastq(char *buf, size_t bsize)
{
    int i, at, newlines;
    for(i = 1, at = -1, newlines = 0; i<bsize; i++)
    {
        if ((buf[i-1] == '\n' ||
             (i>1 && buf[i-1] == '\n' && buf[i-2] == '\r') ||
             (i>1 && buf[i-1] == '\r' && buf[i-2] == '\n'))
             && buf[i] == '@') {
            // it's safe to reset every time a line starts with '@'
            // because if '@' is the first character of the phred string
            // then the next line will be the beginning of the record
            at = i;
            newlines = 0;
        }
        if (buf[i] == '\n')
            newlines++;
        if (at >= 0 &&
            (buf[i-1] == '\n' ||
             (i>1 && buf[i-1] == '\n' && buf[i-2] == '\r') ||
             (i>1 && buf[i-1] == '\r' && buf[i-2] == '\n'))
             && buf[i] == '+')
            // assume that identifier is following (checked sporadically in kvarq.fastq.Fastq)
            break;
    }

    return (at >= 0 && i < bsize) ? at : -1;
}


void scan_filepart(struct scanargs *args)
{
    FILE *fd;
    char *buf, *ptr, *rstart, *rnext, *startread, *plus, *startlongest, *startscore, *seq, *qtr;
    size_t bl;
    long pos, left, recordi;
    int i, j, e, lines, seqi, hits, seql, rl;
    ll_item *tail;
    long buf_recs, buf_tooshort;

    // CAREFUL : opening the file in mode "r" will result in '\r' being discarded
    //           from the bytes read into buffer !
    fd = fopen(args->fname,"rb");
    if (fd == NULL)
    {
        exception = PyExc_MemoryError;
        strncpy(errstr, "cannot open file", ERRSTR_LENGTH);
        return;
    }

    buf = (char *) malloc(SCANBUFSIZE);
    if (buf == NULL)
    {
        exception = PyExc_MemoryError;
        strncpy(errstr, "cannot allocate memory for scanning", ERRSTR_LENGTH);
        fclose(fd);
        return;
    }

    hits = 0;

    tail = args->root; // set to NULL if add_ll can't allocate memory
    recordi = -1;
    // read file chunk buf by buf
    for(pos=args->start,left=args->length; left>0 && !stop;)
    {
        //TODO? mmap instead of read : efficiency vs compatability
        fseek(fd, pos, SEEK_SET);
        bl = fread((void *) buf, 1, MIN(left,SCANBUFSIZE), fd);
        // rnext[0] == '@' (1st byte of next record)
        rnext = buf;

        // DBG("read %li bytes (thread %li)", bl, thread_self());

        // loop over reads in buf
        buf_recs = buf_tooshort = 0;
        while(rnext-buf<bl)
        {
            rstart = rnext;

            if (domore)
                analyse_record(rstart, bl-(rstart-buf));

            // "parse" record
            for(ptr=rstart,lines=0,rl=-1,startread=NULL,startscore=NULL,plus=NULL;
                    lines<4 && ptr-buf<bl;
                    ptr++)
                if (*ptr == '\n')
                {
                    lines++;
                    if (lines == 1)
                        startread = ptr+1;
                    if (lines == 2)
                        plus = ptr+1;
                    if (lines == 3)
                        startscore = ptr+1;
                }

            // don't process partial records
            if (lines<4)
                break;

            // .fastq file format sanity checks
            if (*rstart != '@') {
                exception = fastq_exception;
                snprintf(errstr, ERRSTR_LENGTH, "record must start with '@' (and not '%c') fpos=%ld [junk %ld+%ld]",
                        *rstart, pos + (rstart-buf), args->start, args->length);
                fclose(fd);
                return;
            }
            if (*plus != '+') {
                exception = fastq_exception;
                snprintf(errstr, ERRSTR_LENGTH, "3rd line of record must start with '+' fpos=%ld [junk %ld+%ld]",
                        pos + (plus-buf), args->start, args->length);
                fclose(fd);
                return;
            }

            buf_recs++;

            rnext = ptr;

            // find longest read with good enough quality
            for(ptr=startscore, qtr=startscore, rl=0;
                    ptr==startscore || *(ptr-1)!='\n';
                    ptr++)
                if (*ptr>=Amin) { // '\n' as well as '\r' are <Amin
                    if (!qtr) qtr = ptr;
                } else {
                    if (qtr) {
                        if ((int) (ptr-qtr) > rl) {
                            rl = (int) (ptr-qtr);
                            startlongest = qtr;
                        }
                        qtr = NULL;
                    }
                }
            add_rl(rl);
            startread += startlongest-startscore;

            // dump_record(rstart, startread, rl);

            /*
            if (rl==0) {
                char tmpbuf[1024], *tmpptr;
                memset((void *) tmpbuf, 0, sizeof(tmpbuf));
                for(tmpptr=startscore; *tmpptr!='\n'; tmpptr++)
                    tmpbuf[tmpptr-startscore] = *tmpptr;
                DBG("thread=%ld start=%ld length=%ld", thread_self(), args->start, args->length);
                DBG("rl=%d pos=%ld rstart=%p buf=%p startread=%p startscore=%p ptr=%p qtr=%p",
                        rl, pos, rstart, buf, startread, startscore, ptr, qtr);

                DBG("score=%s", tmpbuf);
            }
            */

            /*
            DBG("parsed : lines=%i startread=%i+%i rl=%i", lines, (int) (startread-buf), (int) (startlongest-startread), (int) rl);
            startread[rl] = 0;
            DBG("  %s", startread);
            startlongest[rl] = 0;
            DBG("  %s", startlongest);
            */

            recordi += 1;
            // DBG("record %li rl=%i fpos(rstart)=%li (thread %li)",
            //         recordi, rl, pos+(rstart-buf), thread_self());

            if (rl<minreadlength) {
                buf_tooshort++;
                continue;
            }

            // DBG("trying record %li (thread %li)", recordi, thread_self());
            // find sequences
            for(seqi=0; args->seqlist[seqi]!=NULL && tail != NULL; seqi++)
            {
                seq = args->seqlist[seqi];
                seql= args->seqlengths[seqi];

                if (rl>minoverlap && seql>minoverlap)
                {
                    // (tail of) read overlaps beginning of sequence
                    // (rl-i<=seql-1) not to count bordercase here and in "read withing seq"
                    for(i=rl-minoverlap; i>0 && rl-i<=seql-1; i--)
                    {
                        for(j=0,e=0; i+j<rl && e<=maxerrors; j++)
                            if (startread[i+j] != seq[j])
                                e++;
                        if (e > maxerrors)
                            continue;
                        hits++;
                        // DBG("adding read where tail overlaps i=%i", i);
                        tail = add_ll(tail, seqi, pos+(startread-buf), -i, rl-i, rl);
                        if (tail == NULL) break;
                    }

                    // (start of) read overlaps end of sequence
                    for(i=seql-minoverlap; i>0 && seql-i<=rl; i--)
                    {
                        for(j=0,e=0; i+j<seql && e<=maxerrors; j++)
                            if (seq[i+j] != startread[j])
                                e++;
                        if (e > maxerrors)
                            continue;
                        hits++;
                        // DBG("adding read where start overlaps i=%i", i);
                        tail = add_ll(tail, seqi, pos+(startread-buf), i, seql-i, rl);
                        if (tail == NULL) break;
                    }
                }

                if (rl>seql)
                {
                    // sequence within read
                    for(i=0; i<=rl-seql; i++)
                    {
                        for(j=0,e=0; j<seql && e<=maxerrors; j++)
                            if (startread[i+j] != seq[j])
                                e++;
                        if (e > maxerrors)
                            continue;
                        hits++;
                        // DBG("adding sequence within read i=%i ", i);
                        tail = add_ll(tail, seqi, pos+(startread-buf), -i, seql, rl);
                        if (tail == NULL) break;
                    }
                }
                else
                {
                    // read within sequence
                    for(i=0; i<=seql-rl; i++)
                    {
                        for(j=0,e=0; j<rl && e<=maxerrors; j++)
                            if (seq[i+j] != startread[j])
                                e++;
                        if (e > maxerrors)
                            continue;
                        hits++;
                        // DBG("adding read within sequence i=%i ", i);
                        tail = add_ll(tail, seqi, pos+(startread-buf), i, rl, rl);
                    }
                }
            }

            if (tail == NULL)
            {
                free(buf);
                exception = PyExc_MemoryError;
                strncpy(errstr, "cannot allocate memory for results", ERRSTR_LENGTH);
                fclose(fd);
                return;
            }

        } // end : loop over reads in buf
        add_records_parsed(buf_recs);

        // incomplete record at end of chunk : break loop
        if (rnext == buf)
            break;

        pos += rnext-buf;
        left -= rnext-buf;
        parsed++;
        //DBG("parsed %li -> %li (parsed=%li) recs=%i tooshort=%i",
        //        pos-(rstart-buf), pos, parsed, buf_recs, buf_tooshort);
    } // end : read file chunk buf by buf

    free(buf);
    args->hits = hits;
    fclose(fd);
}


void distribute_fileparts(struct scanargs *args)
{
    int i = 0;
    while(args[i].fname != NULL) {

        pthread_mutex_lock(&args_mutex);
        if (args[i].assigned == 0) {
            args[i].assigned++;
            pthread_mutex_unlock(&args_mutex);

            DBG("scanning filepart #%i (%li) : %ld..%ld-1",
                    i, thread_self(), args[i].start, args[i].start + args[i].length);
            scan_filepart(args+i);

            // PyErr_Occurred() segfaults because of pthreads ?
            if (exception)
                return;

        } else {
            pthread_mutex_unlock(&args_mutex);
        }
        i++;
    }
}


static PyObject *
engine_findseqs(PyObject *self, PyObject *args)
{
    FILE *fd;
    const char *fname;
    long fsz, length, fpos;
    PyObject *seqlist_obj, *str, *ret, *hits, *stats;
    int i, j, n, threadi, args_n, args_i, err;
    char **seqlist, recbuf[10240];
    int *seqlengths;
    ll_item *root, *item;
    pthread_t *threads;
    struct scanargs *thread_args;

    if (running != 0) //FIXME threadsafe
    {
        PyErr_SetString(PyExc_RuntimeError, "findseqs() already running!");
        return NULL;
    }
    running++;
    stop = 0;
    sigints = 0;

    // argument parsing

    if (!PyArg_ParseTuple(args, "sO", &fname, &seqlist_obj)) {
        running--;
        return NULL;
    }

    if (!PySequence_Check(seqlist_obj))
    {
        PyErr_SetString(PyExc_TypeError, "seqlist must be list of strings");
        running--;
        return NULL;
    }

    seqlist = (char **) malloc((PySequence_Size(seqlist_obj)+1) * sizeof(char *));
    if (seqlist == NULL) {
        running--;
        return PyErr_NoMemory();
    }
    seqlengths = (int *) malloc(PySequence_Size(seqlist_obj) * sizeof(int *));
    if (seqlengths == NULL)
    {
        free(seqlist);
        running--;
        return PyErr_NoMemory();
    }

    for(i=0; i<PySequence_Size(seqlist_obj); i++) 
    {
        str = PySequence_GetItem(seqlist_obj, i);
        seqlist[i] = PyString_AsString(str);
        if (seqlist[i] == NULL)
        {
            free(seqlengths);
            free(seqlist);
            PyErr_SetString(PyExc_TypeError, "seqlist must be list of strings");
            running--;
            return NULL;
        }
        seqlengths[i] = (int) PyString_Size(str);
    }
    seqlist[i] = NULL;

    // prepare sequence quest

    fd = fopen(fname,"rb");
    if (fd == NULL)
    {
        free(seqlengths);
        free(seqlist);
        PyErr_SetString(PyExc_IOError, "cannot open file");
        running--;
        return NULL;
    }
    fseek(fd, 0L, SEEK_END);
    fsz = ftell(fd);

    root = (ll_item *) malloc(sizeof(ll_item));
    if (root == NULL)
    {
        fclose(fd);
        free(seqlist);
        free(seqlengths);
        running--;
        return PyErr_NoMemory();
    }
    root->fpos = -1L;
    root->next = NULL;

    init_stats(seqlist);

    // start sequence quest

    threads = (pthread_t *) malloc(sizeof(pthread_t) *nthreads);
    args_n = MAX(nthreads, MIN(10 * nthreads, fsz / (1024 * 1024 * 10))); // junks >= 10 MB up to n=2*n(threads)
    thread_args = (struct scanargs *) malloc(sizeof(struct scanargs) *(args_n+1));
    total = fsz / SCANBUFSIZE;
    parsed = 0;

    DBG("nthreads=%d args_n=%d", nthreads, args_n);

    for(args_i=0, fpos=0, length=0; args_i<args_n; args_i++)
    {
        if (args_i == args_n-1)
        {
            // last thread gets remainder of file
            length = fsz-fpos > 0 ? fsz-fpos : 0;
        }
        else
        {
            length = fsz/args_n;
            fseek(fd, fpos+length, SEEK_SET);

            //DBG("next record was at : %li", fpos+length);

            // make sure next thread starts at beginning of new record...
            n = fread((void *) recbuf, 1, sizeof(recbuf), fd);
            i = forward_fastq(recbuf, n);
            // ignore error i==-1 : part will fail when scanning

            length += i;

            //DBG("advanced by : %i", i);
        }

        thread_args[args_i].assigned = 0;
        thread_args[args_i].fname = fname;
        thread_args[args_i].start = fpos;
        thread_args[args_i].length = length;
        thread_args[args_i].seqlist = seqlist;
        thread_args[args_i].seqlengths = seqlengths;
        thread_args[args_i].root = root;

        // DBG("thread_args[%d] %ld+%ld", args_i, fpos, length);

        // advance file pointer
        if (length >0)
            fpos += length;
    }
    thread_args[args_i].fname = NULL;

    Py_BEGIN_ALLOW_THREADS

    for(threadi=0; threadi<nthreads; threadi++)
    {

        err = pthread_create(threads+threadi, NULL, 
                (void *(*)(void *)) distribute_fileparts, 
                thread_args);
        DBG("created thread #%d (%li)", threadi, pthread2long(threads[threadi]));

        if (err != 0)
        {
            stop++;
            for(i=0; i<threadi; i++)
                pthread_join(threads[i], NULL);
            exception = PyExc_RuntimeError;
            strncpy(errstr, "pthread_create failed", ERRSTR_LENGTH);
            threadi = -1;
            break;
        }

        // advance file pointer
        if (length >0)
            fpos += length;
    }

    ret = NULL;

    if (threadi != -1)
        for(threadi=0; threadi<nthreads; threadi++)
            pthread_join(threads[threadi], NULL);

    Py_END_ALLOW_THREADS

    if (threadi != -1)
    {

        if (exception == NULL)
        {
            // convert return value

            for(i=0,item=root->next; item!=NULL; i++)
                item = item->next;

            hits = PyTuple_New(i);
            for(i=0,item=root->next; item!=NULL; i++,item=item->next)
                PyTuple_SetItem(hits, i,
                        PyObject_CallObject(hittuple, 
                            Py_BuildValue("(iliii)", 
                                item->seqi,
                                item->fpos,
                                item->spos,
                                item->length,
                                item->readlength)
                            )
                        );

            ret = Py_BuildValue("{sOsO}",
                    "hits", hits,
                    "stats", engine_stats(self, Py_BuildValue("()")));
        }
    }

    // clean up

    free(threads);
    free(thread_args);
    free_ll(root);
    fclose(fd);
    free(seqlist);
    free(seqlengths);

    if (exception != NULL) {
        PyErr_SetString(exception, errstr);
        exception = NULL;
    }

    running--;
    return ret;
}



static PyObject *
engine_get_config(PyObject *self, PyObject *args)
{
    return Py_BuildValue("{sisisisiscsc}", 
            "maxerrors", maxerrors,
            "minoverlap", minoverlap,
            "minreadlength", minreadlength,
            "nthreads", nthreads,
            "Amin", Amin,
            "Azero", Azero);
}


static PyObject *
engine_config(PyObject *self, PyObject *args, PyObject *kw)
{
    static char *kwl[] = { "maxerrors", "minoverlap", "minreadlength", "nthreads", "Amin", "Azero", NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kw, 
                "|iiiicc", kwl, &maxerrors, &minoverlap, &minreadlength, &nthreads, &Amin, &Azero))
        return NULL;

    Py_RETURN_NONE;
}


static PyObject *
engine_abort(PyObject *self, PyObject *args)
{
    stop++;
    Py_RETURN_NONE;
}


void *do_nothing(void *ptr) {
    int i=0, x, y;
    for(x=0; x<100000; x++)
    for(y=0; y<30000; y++) {
        i-= i%15; // cannot be optimized away easily...
        i+= (i/2 +5) *3;
    }
    return (void *) 0;
}

static PyObject *
engine_loop(PyObject *self, PyObject *args)
{
    int n = 8, i;
    pthread_t *pts;
    int *prs;

    pts = (pthread_t *) malloc(sizeof(pthread_t) *n);
    prs = (int       *) malloc(sizeof(int      ) *n);
    for(i=0; i<n; i++)
        prs[i] = pthread_create(pts+i, NULL, do_nothing, NULL);

    for(i=0; i<n; i++)
        pthread_join(pts[i], NULL);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
engine_test(PyObject *self, PyObject *args)
{
    PyObject *number, *tuple;
    
    number = PyInt_FromLong(1L);
    tuple = PyTuple_New(1);
    printf("number->refcount=%d\n", -1);
    Py_RETURN_NONE;
}


static PyMethodDef methods[] = {
    {"test",  engine_test, METH_VARARGS,
     "Perform some test."},
    {"config",  (PyCFunction)engine_config, METH_VARARGS | METH_KEYWORDS,
     "configure the engine; keywords are:\n\n"
     "maxerrors : maximum number of base mismatches in sequence alignment\n"
     "minoverlap : minimum number of base overlap for beginning/end hits\n"
     "minreadlength : ignore reads shorter than this\n"
     "nthreads : number of threads to use for scanning\n"
     "Amin : nucleotides with quality ASCII value lower than this are discarded\n"
     "Azero : ASCII value that corresponds to Q=0 (depends on FastQ format)\n"},
    {"get_config", engine_get_config, METH_VARARGS,
     "get the current config as dictionary"},
    {"abort",  engine_abort, METH_VARARGS,
     "Stop the current scanning."},
    {"loop",  engine_loop, METH_VARARGS,
     "Do nothing, but real hard."},
    {"findseqs", engine_findseqs, METH_VARARGS,
     "finds occurences of base sequences in fastq files\n\n"
     "fname : filename of fastq file\n"
     "sequences : list of sequences to look for\n\n"
     "returns a dictionary 'hits' and 'stats':\n\n"
     "'hits' is a named tuple of (seq_nr, file_pos, seq_pos, length, readlength)\n"
     "file_pos describes beginning of read, seq position places the\n"
     "beginning of the read relative to the beginning of the sequence\n"
     "(<0 if read overlaps only with beginning of sequence or read contains\n"
     "whole sequence; >0 if read overlaps only with end of sequence or read\n"
     "is contained within sequence), length gives the number of overlapping basepairs\n\n"
     "'stats' is the same dict as returned by a call to stats()\n"},
    {"stop",  engine_stop, METH_VARARGS,
     "stop the scanning process"},
    {"stats", engine_stats, METH_VARARGS,
     "returns a dict containing information about scanning process\n"
     "  - 'readlengths' : tuple of number of occurences when accessed by read length\n"
     "  - 'progress' : current progress of findseqs() ranging 0..1\n"
     "  - 'sigints' : how many <CTRL-C> were caught since beginning of scan\n"
     "  - 'nseqbasehits' : sum(hit_length), indexed by sequence as given to findseqs()\n"
     "  - 'nseqhits' : number of hits, indexed by sequence as given to findseqs()\n"
     "  - 'records_parsed' : total number of records parsed\n"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


// if called initX : "module does not defined init function (initc)"
// if called initd and setup.py calls extensions 'd' : "dynamic module not initialized properly"
PyMODINIT_FUNC
initengine(void)
{
    PyObject *collections, *namedtuple, *fastq;

    (void) Py_InitModule("engine", methods);

    collections = PyImport_ImportModule("collections");
    namedtuple = PyObject_GetAttrString(collections, "namedtuple");

    hittuple = PyObject_CallObject(namedtuple,
            Py_BuildValue("(sO)", "Hit", 
                    Py_BuildValue("(sssss)",
                            "seq_nr",
                            "file_pos",
                            "seq_pos",
                            "length",
                            "readlength")));
    engine_mod = PyImport_ImportModule("kvarq.engine");
    PyObject_SetAttrString(engine_mod, "Hit", hittuple);

    fastq = PyImport_ImportModule("kvarq.fastq");
    fastq_exception = PyObject_GetAttrString(fastq, "FastqFileFormatException");
    exception = NULL;

#ifdef _WIN32 // 32/64 bit windows
    SetConsoleCtrlHandler((PHANDLER_ROUTINE) CtrlHandler, TRUE);
#else // POSIX
    signal(SIGINT, sigint_cb);
#endif
}

