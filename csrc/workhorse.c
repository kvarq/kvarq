#include <Python.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <stdarg.h>

#include "gz/miniz.c"

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


/* data structures {{{1 */

typedef struct ll_item_ {
    int seqi;
    long fpos;
    int spos;
    int length;
    int readlength;
    struct ll_item_ *next;
} ll_item;

struct fastq_file {
    const char **fnames; // NULL terminated
    int fname_i;
    const char *buf; // partial record from last read
    size_t buf_size;
    FILE *fd;
    size_t size;
    size_t fpos; // within inflated data
    size_t ftell0; // sum of filesizes of files already read
    int eof;

    int compressed;
    mz_stream mzs;
    char *inbuf;
    long remaining;
};

struct scanargs {
    struct fastq_file *fastq;
    char **seqlist;
    int *seqlengths;
    ll_item *root;
    PyObject *pyhitseqs;
};


/* globals {{{1 */

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
pthread_mutex_t fastq_read_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t records_parsed_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t interpreter_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t profile_mutex = PTHREAD_MUTEX_INITIALIZER;
// python objects for interfacing etc
PyObject *engine_mod, *hittuple;
// python pthreads exception handling
PyObject *exception, *fastq_exception;
// python object for logging output
PyObject *lo_log;
long LOG_DEBUG, LOG_INFO, LOG_WARNING, LOG_ERROR, LOG_FATAL;
#define LOG_BUF_SIZE 1024
#define ERRSTR_LENGTH 1024
char errstr[ERRSTR_LENGTH];

// stats
size_t fastq_size_estimated, fastq_parsed;
int sigints, nseqs;
long *seqbasehits, *seqhits, records_parsed;
long nrecords;
long nN, nG, nA, nC, nT, nX;
#define MAX_READLENGTH 1024
long rls_longest;
long rls_buf[MAX_READLENGTH];
#define AMIN_STEPS 5
long **all_rls_buf;


/* signal handler {{{1 */

#ifdef _WIN32 // 32/64 bit windows
BOOL CtrlHandler(DWORD fdwCtrlType)
{
    // DBG("CtrlHandler called");
    if (fdwCtrlType == CTRL_C_EVENT)
    {
	// DBG("CtrlHandler detected CTRL_C_EVENT");
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


/* utility functions {{{1 */

// use globals LOG_{DEBUG|INFO|WARNING|ERROR|FATAL} as first argument

void lo_log_msg(int level, char *fmt, ...)
{
    char *log_header = "[kvarq.engine] ", buf[LOG_BUF_SIZE];
    va_list args;

    sprintf(buf, "%s", log_header);
    va_start(args, fmt);
    vsnprintf(buf + strlen(buf), LOG_BUF_SIZE - strlen(buf), fmt, args);
    va_end(args);

    pthread_mutex_lock(&interpreter_mutex);
    PyObject_CallObject(lo_log, Py_BuildValue("(is)", level, buf));
    pthread_mutex_unlock(&interpreter_mutex);
}

void dump_record(char *startrecord, char *startread, int rl)
{
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


/* thread profiling {{{1 */

#if 0

struct profile_struct {
    char *name;
    long threadid;
    clock_t clocks, running;
};

#define PROFILES_SIZE_STEP 100
int profiles_size = 0, profiles_n = 0;
struct profile_struct *profiles = NULL;

void profile_start(char *name)
{
    int i;
    long threadid = thread_self();

    pthread_mutex_lock(&profile_mutex);

    for(i = 0; i < profiles_n; i++)
	if (strcmp(profiles[i].name, name) == 0 &&
		profiles[i].threadid == threadid)
	    break;

    if (i == profiles_n)
    {
	if (profiles_n == profiles_size)
	{
	    profiles_size += PROFILES_SIZE_STEP;
	    profiles = realloc(profiles, sizeof(struct profile_struct) * profiles_size);
	}
	profiles[i].name = name;
	profiles[i].threadid = threadid;
	profiles[i].clocks = 0;
	profiles[i].running = 0;

	profiles_n++;
    }

    profiles[i].running = clock();

    pthread_mutex_unlock(&profile_mutex);
}

void profile_stop(char *name)
{
    int i;
    long threadid = thread_self();

    pthread_mutex_lock(&profile_mutex);

    for(i = 0; i < profiles_n; i++)
	if (strcmp(profiles[i].name, name) == 0 &&
		profiles[i].threadid == threadid)
	    break;

    profiles[i].clocks += clock() - profiles[i].running;
    profiles[i].running = 0;

    pthread_mutex_unlock(&profile_mutex);
}

void profile_dump()
{
    int i;

    fputc('\n', stderr);
    for(i = 0; i < profiles_n; i++)
	fprintf(stderr, "[%ld] %12s clocks=%10ld running=%ld\n",
		profiles[i].threadid, profiles[i].name,
		profiles[i].clocks, profiles[i].running);
    fputc('\n', stderr);
}

#else

int profiles_n = 0;

void profile_start(char *name) {}
void profile_stop(char *name) {}
void profile_dump() {}

#endif


/* stats functions {{{1 */

/* init {{{2 */

void init_stats(char **seqlist)
{
    int i;

    //FIXME : allocated structures are never freed

    memset((void *) rls_buf, 0, sizeof(rls_buf));
    rls_longest = -1;

    nrecords = nA = nG = nC = nT = nN = nX = 0;

    if (all_rls_buf) {
	for(i=0; i<2*AMIN_STEPS; i++)
	    free(all_rls_buf[i]);
	free(all_rls_buf);
    }
    all_rls_buf = (long **) malloc(sizeof(long *) * AMIN_STEPS*2);
    for(i=0; i<2*AMIN_STEPS; i++) {
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

/* add infos {{{2 */

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
    for(i=0; i<2*AMIN_STEPS; i++)
	for(++ptr, qtr=ptr; *ptr!='\n'; ptr++)
	    if (*ptr<(Amin+(i<AMIN_STEPS?-i-1:i-AMIN_STEPS+1))) {
		if (qtr)
		    all_rls_buf[i][qtr-ptr>=MAX_READLENGTH ? MAX_READLENGTH-1 : ptr-qtr]++;
		qtr = NULL;
	    } else if (!qtr) qtr=ptr;
}

void add_records_parsed(long n) {
    pthread_mutex_lock(&records_parsed_mutex);
    // DBG("records_parsed : %ld -> %ld", records_parsed, records_parsed + n);
    records_parsed += n;
    pthread_mutex_unlock(&records_parsed_mutex);
}

void add_rl(long rl)
{
    pthread_mutex_lock(&rl_mutex);
    if (rl>=0 && rl<MAX_READLENGTH)
	rls_buf[rl]++;
    if (rl > rls_longest)
	rls_longest = rl;
    pthread_mutex_unlock(&rl_mutex);
}

/* findseqs {{{1 */

/* adding hits {{{2 */

ll_item *add_hit(ll_item *item, int seqi, long fpos, int spos, int length, int readlength, PyObject *pyhitseqs, char *hitseq)
{
    PyObject *pyhitseq;
    pthread_mutex_lock(&ll_mutex);

    //DBG("adding item seqi=%d fpos=%li spos=%i length=%i (thread %li)",
    //        seqi, fpos, spos, length, thread_self());

    while(item->next != NULL)
	item = item->next;
    item->next = (ll_item *) malloc(sizeof(ll_item));

    if (item->next == NULL)
    {
	pthread_mutex_unlock(&ll_mutex);
	return (ll_item *) PyErr_NoMemory();
    }

    item->next->seqi = seqi;
    item->next->fpos = fpos;
    item->next->spos = spos;
    item->next->length = length;
    item->next->readlength = readlength;
    item->next->next = NULL;
    // DBG("added %ld at %ld (thread %ld)", which, pos, thread_self());

    seqbasehits[seqi] += length;
    seqhits[seqi]++;

    pyhitseq = PyString_FromStringAndSize(hitseq, length);
    PyList_Append(pyhitseqs, pyhitseq);
    Py_XDECREF(pyhitseq);

    pthread_mutex_unlock(&ll_mutex);
    return item->next;
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

/* parse .gz {{{2 */

#define GZ_DEFLATED     8
#define GZ_ASCII_FLAG   0x01 /* bit 0 set: file probably ascii text */
#define GZ_CONTINUATION 0x02 /* bit 1 set: continuation of multi-part gzip file */
#define GZ_EXTRA_FIELD  0x04 /* bit 2 set: extra field present */
#define GZ_ORIG_NAME    0x08 /* bit 3 set: original file name present */
#define GZ_COMMENT      0x10 /* bit 4 set: file comment present */
#define GZ_ENCRYPTED    0x20 /* bit 5 set: file is encrypted */
#define GZ_RESERVED     0xC0 /* bit 6,7:   reserved */

/**
 * advances file position past gz header
 *
 * checks for a valid gz header (starting with magix bytes 0x1F 0x8B)
 * in file stream and advances file pointer past the header; no data
 * is returned but the correctness of the gz format (without encryption,
 * continuation, or reserved features) is asserted
 *
 * @param fd FILE pointer
 * @param dist how many random bytes can be accepted before the header
 *
 * @return NULL if successful, otherwise descriptive error string
 */

const char *skip_gz_header(FILE *fd, int dist)
{
    // (based on gzip 1.2.4 source code)
    // 2 bytes : magic 1F 8B
    // 1 byte  : method -- must be DEFLATED
    // 1 byte  : flags -- does not support CONTINUATION, ENCRYPTED, RESERVED
    // 4 bytes : timestamp
    // 1 byte  : extra flags
    // 1 byte  : os type
    // optional : original name
    // optional : comment

    int state, c, flags, i, n, y;

    for(state=0, c=fgetc(fd), y=0; state != 2 && y<=dist && c != -1; c=fgetc(fd)) {
        if (c == 0x1F && state == 0)
            state++;
        else if (c == 0x8B && state ==1)
            state++;
        else {
            state=0;
            y++;
        }
    }

    if (state != 2)
        return "magic bytes not found";
    //if (y) fprintf(stderr, "ignored %d<%d prior to header\n", y, dist);

    if (c != GZ_DEFLATED) {
        return "expected method==DEFLATED";
    }
    flags = fgetc(fd);
    if (flags & (GZ_CONTINUATION | GZ_ENCRYPTED | GZ_RESERVED)) {
        return "unsupported flags (CONTINUATION or ENCRYPTED or RESERVED)";
    }
    for(i=0; i<4 + 2; i++)
        // ignore stamp, extra flags, os type
        (void) fgetc(fd);
    if (flags & GZ_EXTRA_FIELD) {
        n = fgetc(fd);
        n|= fgetc(fd) << 8;
        // fprintf(stderr, "ingoring extra field length %d\n", n);
        while(n--) (void) fgetc(fd);
    }
    if (flags & GZ_ORIG_NAME) {
        // fprintf(stderr, "ignoring original name\n");
        do {
            c = fgetc(fd);
        } while (c>0);
    }
    if (flags & GZ_COMMENT) {
        // fprintf(stderr, "ignoring comment\n");
        do {
            c = fgetc(fd);
        } while (c>0);
    }

    return NULL; // success
}


/* read from .fastq files {{{2 */

/**
 * opens next file in sequence of fastq->fnames
 *
 * this next file can be a .fastq or .fastq.gz file regardless of
 * whether the last file type (although time estimates will be totally
 * wrong if filetypes are mixed)
 *
 * @param fastq fastq_structure where fastq_i points to file that
 *     should be opened next
 * @return 0 in case of success, -1 in case of error (globals exception
 *     and errstr are set accordingly)
 */

int fastq_open_next(struct fastq_file *fastq)
{
    const char *fname, *ret;

    if (fastq->fname_i > 0) {
	// close open file & add bytes already read
	fastq->ftell0 += ftell(fastq->fd);
	fclose(fastq->fd);
    }

    fname = fastq->fnames[fastq->fname_i];
    fastq->fname_i++;

    // CAREFUL : opening the file in mode "r" will result in '\r' being discarded
    //           from the bytes read into buffer !
    fastq->fd = fopen(fname,"rb");
    if (fastq->fd == NULL)
    {
	exception = PyExc_IOError;
	snprintf(errstr, ERRSTR_LENGTH, "cannot open file");
	return -1;
    }

    if (strcmp(fname + strlen(fname) - 3, ".gz") == 0)
    {
	// initialize datastructure for inflating
	fastq->compressed = 1;

	if (mz_inflateInit2(&fastq->mzs, -MZ_DEFAULT_WINDOW_BITS) != MZ_OK)
	{
	    fclose(fastq->fd);
	    exception = PyExc_RuntimeError;
	    snprintf(errstr, ERRSTR_LENGTH, "cannot mz_inflateInit()");
	    return -1;
	}

	if (fastq->inbuf == NULL)
	    fastq->inbuf = malloc(SCANBUFSIZE);

	if (fastq->inbuf == NULL)
	{
	    fclose(fastq->fd);
	    exception = PyExc_MemoryError;
	    snprintf(errstr, ERRSTR_LENGTH, "cannot allocate inbuf");
	    return -1;
	}

	fastq->mzs.next_in = (const unsigned char *) fastq->inbuf;
	fseek(fastq->fd, 0, SEEK_END);
	fastq->remaining = ftell(fastq->fd);
	fseek(fastq->fd, 0, SEEK_SET);
	ret = skip_gz_header(fastq->fd, 0);
	if (ret != NULL)
	{
	    fclose(fastq->fd);
	    exception = PyExc_IOError;
	    snprintf(errstr, ERRSTR_LENGTH, "no valid gzip header found "
		    "at beginning of file : %s", ret);
	    return -1;
	}
	fastq->remaining -= ftell(fastq->fd);

	// random guess
	fastq_size_estimated *= 3;
    }

    return 0;
}

/**
 * opens a .fastq file for further access via fastq_read
 *
 * also initializes globals fastq_size_estimated and fastq_parsed
 *
 * @param fnames NULL terminated array of paths of the .fastq files
 * @return pointer to fastq file object or NULL in case of error
 *         (PyErr_SetString called with appropriate arguments)
 */

struct fastq_file *fastq_open(const char **fnames)
{
    struct fastq_file *fastq;
    int i;
    FILE *fd;

    fastq = (struct fastq_file *) malloc(sizeof(struct fastq_file));
    if (fastq == NULL)
    {
	exception = PyExc_MemoryError;
	snprintf(errstr, ERRSTR_LENGTH, "cannot allocate struct fastq_file");
	return NULL;
    }

    memset(fastq, 0, sizeof(struct fastq_file));
    fastq->fnames = fnames;

    // determine overall file size
    for(i = 0; fastq->fnames[i]; i++)
    {
	fd = fopen(fastq->fnames[i], "rb");
	if (fd == NULL)
	{
	    free(fastq);
	    exception = PyExc_IOError;
	    snprintf(errstr, ERRSTR_LENGTH,
		    "cannot open file '%s' for getting filesize", fastq->fnames[i]);
	    return NULL;
	}
	fseek(fd, 0, SEEK_END);
	fastq->size += ftell(fd);
	fclose(fd);
    }

    // initialize globals fastq_parsed, fastq_size_estimated
    fastq_parsed = 0;
    fastq_size_estimated = fastq->size;

    if (fastq_open_next(fastq) != 0)
    {
	free(fastq);
	return NULL;
    }

    return fastq;
}

/**
 * get length of last partial fastq record in buffer
 *
 * @param buf data read from .fastq file
 * @param n length of data in buf
 * @return length of last partial record in buf; -1 if no partial buf found at end
 */

int fastq_rewind(char *buf, int n)
{

    // we scan for /^\+/ (line3) and then accept the next /^@/ (line1) as the
    // beginning of a record note that /^\+/ could also be a PHRED score
    // therefore accept a repeated /^\+/ as the "real" line3

    int state = 0; // 1=found line3
    int i = 1;

    while (i + 1 < n)
    {
	if (buf[n - i] == '+' && 
		(buf[n - i - 1] == '\n' || buf[n - i - 1] == '\r'))
	    state = 1; // state might have been 0 OR 1 (see note above)
	else if (state == 1 && buf[n - i] == '@' &&
		(buf[n - i - 1] == '\n' || buf[n - i - 1] == '\r'))
	    return i;
	i++;
    }

    return -1;
}

/**
 * reads content from a .fastq file -- multithread safe
 *
 * also updates globals fastq_size_estimated and fastq_parsed
 *
 * @param fastq pointer as returned by fastq_open
 * @param buf where the read data is aved
 * @param buf_size maximum number of bytes that can be saved in buf
 * @param fposp pointer to file position of first byte in buffer
 *
 * @return number of bytes read; this is normally smaller than buf_size
 *         and the last character is guaranteed to be the end of a
 *         record (unless the file ends with a partial record).
 *         returns 0 if EOF is encountered and -1 in case of error
 *         (with exception/errstr globals accordingly set)
 */

size_t fastq_read(struct fastq_file *fastq, char *buf, size_t buf_size, size_t *fposp)
{
    size_t n, leftovers, m;
    long pos;
    int status;
    unsigned int avail_in, avail_out;
    const char *msg;

    profile_start("fastq_read");
    pthread_mutex_lock(&fastq_read_mutex);

    // eof? open next file if available
    if (fastq->eof != 0)
	if (fastq->fnames[fastq->fname_i] != NULL)
	    if (fastq_open_next(fastq) != 0)
		return -1;

    // copy partial record from last read
    if (fastq->buf_size > 0)
    {
	if (fastq->buf_size > buf_size) {
	    exception = PyExc_RuntimeError;
	    strncpy(errstr, "buf_size < fastq->buf_size !", ERRSTR_LENGTH);

	    pthread_mutex_unlock(&fastq_read_mutex);
	    profile_stop("fastq_read");
	    return -1;
	}

	leftovers = fastq->buf_size;
	memcpy(buf, fastq->buf, leftovers);
	free((void *) fastq->buf);
	fastq->buf_size = 0;

    } else {
	leftovers = 0;
    }

    // first byte of buffer within (inflated) file
    *fposp = fastq->fpos - leftovers;
    n = 0;
    fastq->eof = 0;

    if (fastq->compressed)
    {
	// read compressed data
	fastq->mzs.next_out = (unsigned char *) (buf + leftovers);
	fastq->mzs.avail_out = buf_size - leftovers;
	do
	{
	    // fill up (deflated) inbuf if empty
	    if (fastq->mzs.avail_in == 0)
	    {
		m = MIN(SCANBUFSIZE, fastq->remaining);
		if (fread(fastq->inbuf, 1, m, fastq->fd) != m)
		{
		    exception = PyExc_IOError;
		    strncpy(errstr, "could not read enough bytes from .fastq.gz", ERRSTR_LENGTH);
		    if (ferror(fastq->fd) != 0)
			strcat(errstr, " : I/O error");
		    if (feof(fastq->fd) != 0)
			strcat(errstr, " : premature EOF");

		    pthread_mutex_unlock(&fastq_read_mutex);
		    profile_stop("fastq_read");
		    return -1;
		}
		fastq->mzs.next_in = (const unsigned char *) fastq->inbuf;
		fastq->mzs.avail_in = m;
		fastq->remaining -= m;
	    }

	    // inflate inbuf -> buf
	    avail_in = fastq->mzs.avail_in;
	    avail_out = fastq->mzs.avail_out;
	    //DBG("will inflate avail_in=%d avail_out=%d", avail_in, avail_out);
	    status = mz_inflate(&fastq->mzs, Z_SYNC_FLUSH);

	    if ((status != MZ_OK) && (status != MZ_STREAM_END))
	    {
		exception = PyExc_IOError;
		snprintf(errstr, ERRSTR_LENGTH,
			"error while inflating compressed data : status=%d"
			" fpos=%ld ftell=%ld+%ld avail_in=%d",
			status, fastq->fpos, fastq->ftell0, ftell(fastq->fd),
			fastq->mzs.avail_in);
		// strncpy(errstr, "error while inflating compressed data", ERRSTR_LENGTH);

		pthread_mutex_unlock(&fastq_read_mutex);
		profile_stop("fastq_read");
		return -1;
	    }

	    // update bytes read in deflated stream
	    n += avail_out - fastq->mzs.avail_out;

	    /*
	    DBG("fastq_read [%li] : inflated %d -> %d bytes",
		    thread_self(),
		    avail_in - fastq->mzs.avail_in, avail_out - fastq->mzs.avail_out);
	    */

	    // some versions of gzip create files with many successive deflated streams
	    if (status == MZ_STREAM_END &&
		    fastq->remaining + fastq->mzs.avail_in > 10)
	    {
		// rewind file
		fseek(fastq->fd, -((long) fastq->mzs.avail_in), SEEK_CUR);
		fastq->remaining += fastq->mzs.avail_in;
		// read header -- 2-7 bytes between streams !?
		pos = ftell(fastq->fd);
		msg = skip_gz_header(fastq->fd, 10);
		if (msg != NULL)
		{
		    lo_log_msg(LOG_ERROR, "cannot read next deflated stream "
			    "in compressed file : %s", msg);
		    fastq->remaining = 0;
		}
		else
		{
		    fastq->remaining -= ftell(fastq->fd) - pos;
		    // reset stream
		    mz_inflateEnd(&fastq->mzs);
		    mz_inflateInit2(&fastq->mzs, -MZ_DEFAULT_WINDOW_BITS);
		    // reset input
		    fastq->mzs.avail_in = 0;
		}
	    }

	} while(fastq->remaining + fastq->mzs.avail_in > 10 &&
		fastq->mzs.avail_out > 0);

	if (fastq->remaining + fastq->mzs.avail_in < 10)
	    fastq->eof = 1;

	//DBG("status=%d mzs.avail_out=%ld", status, fastq->mzs.avail_out);

	// update guess
	fastq_size_estimated = (size_t) (
		(float) fastq->size * (fastq->fpos + n) / (fastq->ftell0 + ftell(fastq->fd)));
    }
    else
    {
	// read uncompressed data
	//TODO? mmap instead of read : efficiency vs compatability
	n += fread((void *) (buf + leftovers), 1, buf_size - leftovers, fastq->fd);

	if (ferror(fastq->fd) != 0) {
	    exception = PyExc_IOError;
	    strncpy(errstr, "error while reading from file in fastq_read", ERRSTR_LENGTH);

	    pthread_mutex_unlock(&fastq_read_mutex);
	    profile_stop("fastq_read");
	    return -1;
	}

	if (feof(fastq->fd) != 0)
	    fastq->eof = 1;
    }

    if (n == 0) {
	// DBG("fastq_read [%li] : reached end of file", thread_self());
	pthread_mutex_unlock(&fastq_read_mutex);
	profile_stop("fastq_read");
	return leftovers + n;
    }

    // DBG("updating filepos %ld -> %ld", fastq->fpos, fastq->fpos + n);
    fastq->fpos += n;
    fastq_parsed += n;

    if (fastq->eof == 0)
    {
	// save partial record at end to bufer
	fastq->buf_size = fastq_rewind(buf, leftovers + n);

	if (fastq->buf_size == -1) {
	    exception = PyExc_RuntimeError;
	    snprintf(errstr, ERRSTR_LENGTH,
		    "could find beginning of record; read %ld bytes up to %ld",
		    n, ftell(fastq->fd));

	    pthread_mutex_unlock(&fastq_read_mutex);
	    profile_stop("fastq_read");
	    return -1;
	}
	if (fastq->buf_size > 0) {
	    fastq->buf = (char *) malloc(fastq->buf_size);
	    if (fastq->buf == NULL)
	    {
		exception = PyExc_MemoryError;
		PyErr_SetString(PyExc_MemoryError, "cannot allocate new fastq->buf");
		profile_stop("fastq_read");
		return -1;
	    }
	    memcpy((void *) fastq->buf, buf + leftovers + n - fastq->buf_size,
		    fastq->buf_size);
	}
    }

    /*
    DBG("fastq_read [%li] : read %ld=(%ld+%ld+%ld) bytes -> fpos=%ld",
	    thread_self(), leftovers + n + fastq->buf_size,
	    leftovers, n, fastq->buf_size,
	    fastq->fpos);
    */

    pthread_mutex_unlock(&fastq_read_mutex);
    profile_stop("fastq_read");
    return leftovers + n - fastq->buf_size;
}

void fastq_close(struct fastq_file *fastq)
{
    fclose(fastq->fd);
    if (fastq->buf_size)
	free((void *) fastq->buf);
    if (fastq->inbuf != NULL)
	free((void *) fastq->inbuf);
    free(fastq);
}

/* scan_filepart {{{2 */

/**
 * scans .fastq file one BUFSIZE at a time
 *
 * sets globals exception, errstr if exception occured
 */

void scan_filepart(struct scanargs *args)
{
    char *buf, *ptr, *rstart, *rnext, *startread, *plus, *startlongest, *startscore, *seq, *qtr;
    size_t bl;
    size_t fpos;
    long recordi;
    int i, j, e, lines, seqi, seql, rl;
    ll_item *tail;
    long buf_recs, buf_tooshort;

    buf = (char *) malloc(SCANBUFSIZE);
    if (buf == NULL)
    {
	exception = PyExc_MemoryError;
	strncpy(errstr, "cannot allocate memory for scanning", ERRSTR_LENGTH);
	return;
    }

    tail = args->root; // set to NULL if add_ll can't allocate memory
    recordi = -1;

    // read file buf by buf {{{3
    while((bl = fastq_read(args->fastq, buf, SCANBUFSIZE, &fpos)) > 0
	    && exception == NULL && stop == 0)
    {
	profile_start("scan buf");

	// rnext[0] == '@' (1st byte of next record)
	rnext = buf;

	// DBG("read %li bytes (thread %li)", bl, thread_self());

	// loop over reads in buf {{{4
	buf_recs = buf_tooshort = 0;
	while(rnext - buf < bl)
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
		snprintf(errstr, ERRSTR_LENGTH, "record must start with '@' (and not '%c') "
			"fpos=%ld", *rstart, fpos + (rstart-buf));
		return;
	    }
	    if (*plus != '+') {
		exception = fastq_exception;
		snprintf(errstr, ERRSTR_LENGTH, "3rd line of record must start with '+' fpos=%ld",
			fpos + (plus-buf));
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
	       DBG("rl=%d fpos=%ld rstart=%p buf=%p startread=%p startscore=%p ptr=%p qtr=%p",
	       rl, fpos, rstart, buf, startread, startscore, ptr, qtr);

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
	    //    recordi, rl, fpos+(rstart-buf), thread_self());

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
			// DBG("adding read where tail overlaps i=%i", i);
			tail = add_hit(tail, seqi, fpos+(startread-buf), -i, rl-i, rl,
				args->pyhitseqs, startread + i);
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
			// DBG("adding read where start overlaps i=%i", i);
			tail = add_hit(tail, seqi, fpos+(startread-buf), i, seql-i, rl,
				args->pyhitseqs, startread);
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
			// DBG("adding sequence within read i=%i ", i);
			tail = add_hit(tail, seqi, fpos+(startread-buf), -i, seql, rl,
				args->pyhitseqs, startread + i);
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
			// DBG("adding read within sequence i=%i ", i);
			tail = add_hit(tail, seqi, fpos+(startread-buf), i, rl, rl,
				args->pyhitseqs, startread);
		    }
		}
	    }

	    if (tail == NULL)
	    {
		free(buf);
		exception = PyExc_MemoryError;
		strncpy(errstr, "cannot allocate memory for results", ERRSTR_LENGTH);
		return;
	    }

	}
	// end : loop over reads in buf }}}4
	add_records_parsed(buf_recs);

	//DBG("parsed %li -> %li (parsed=%li) recs=%i tooshort=%i",
	//        pos-(rstart-buf), pos, parsed, buf_recs, buf_tooshort);
	profile_stop("scan buf");
    }
    // end : read file buf by buf }}}3

    // exception, errstr set if bl == -1
    free(buf);
}


/* module functions {{{1 */

/* engine.stats {{{2 */

    static PyObject *
engine_stats(PyObject *self, PyObject *args)
{
    int i;
    PyObject *rls, *sbhs, *shs;
    float progress;

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

    progress = 0;
    if (fastq_size_estimated > 0)
	progress = ((float) MIN(fastq_parsed, fastq_size_estimated)) / fastq_size_estimated;

    return Py_BuildValue("{sOsfsOsOsOsOsOsO}",
	    "readlengths", rls,
	    "progress", progress,
	    "nseqbasehits", sbhs,
	    "nseqhits", shs,
	    "parsed", PyInt_FromLong(fastq_parsed),
	    "total", PyInt_FromLong(fastq_size_estimated),
	    "sigints", PyInt_FromLong(sigints),
	    "records_parsed", PyInt_FromLong(records_parsed)
	    );
}

/* engine.findseqs {{{2 */

    static PyObject *
engine_findseqs(PyObject *self, PyObject *findseqs_args)
{
    const char **fnames;
    PyObject *fname_obj, *seqlist_obj, *str, *ret, *hits, *pystats;
    int i, threadi, err;
    ll_item *item;
    pthread_t *threads;
    struct scanargs args;

    if (running != 0) //FIXME threadsafe
    {
	PyErr_SetString(PyExc_RuntimeError, "findseqs() already running!");
	return NULL;
    }
    running++;
    stop = 0;
    sigints = 0;

    // argument parsing {{{3

    if (!PyArg_ParseTuple(findseqs_args, "OO", &fname_obj, &seqlist_obj)) {
	running--;
	return NULL;
    }

    if (PyString_Check(fname_obj)) {
	fnames = (const char **) malloc(sizeof(char *) * 2);
	if (fnames == NULL) {
	    running--;
	    return PyErr_NoMemory();
	}
	fnames[0] = PyString_AsString(fname_obj);
	fnames[1] = NULL;

    } else if (PySequence_Check(fname_obj)) {
	fnames = (const char **) malloc(sizeof(char *) *
		(PySequence_Size(fname_obj) + 1));
	if (fnames == NULL) {
	    running--;
	    return PyErr_NoMemory();
	}
	for(i = 0; i < PySequence_Size(fname_obj); i++)
	    fnames[i] = PyString_AsString(PySequence_GetItem(fname_obj, i));
	fnames[i] = NULL;

    } else {
	PyErr_SetString(PyExc_TypeError, "fname must be [sequence of] string[s]");
	running--;
	return NULL;
    }

    if (!PySequence_Check(seqlist_obj))
    {
	PyErr_SetString(PyExc_TypeError, "seqlist must be sequence of strings");
	free(fnames);
	running--;
	return NULL;
    }

    args.seqlist = (char **) malloc((PySequence_Size(seqlist_obj)+1) * sizeof(char *));
    if (args.seqlist == NULL) {
	free(fnames);
	running--;
	return PyErr_NoMemory();
    }
    args.seqlengths = (int *) malloc(PySequence_Size(seqlist_obj) * sizeof(int *));
    if (args.seqlengths == NULL)
    {
	free(args.seqlist);
	free(fnames);
	running--;
	return PyErr_NoMemory();
    }

    for(i=0; i<PySequence_Size(seqlist_obj); i++) 
    {
	str = PySequence_GetItem(seqlist_obj, i);
	args.seqlist[i] = PyString_AsString(str);
	if (args.seqlist[i] == NULL)
	{
	    free(args.seqlengths);
	    free(args.seqlist);
	    PyErr_SetString(PyExc_TypeError, "seqlist must be list of strings");
	    free(fnames);
	    running--;
	    return NULL;
	}
	args.seqlengths[i] = (int) PyString_Size(str);
    }
    args.seqlist[i] = NULL;


    // prepare sequence quest {{{3

    args.fastq = fastq_open(fnames);
    if (args.fastq == NULL)
    {
	PyErr_SetString(exception, errstr);
	free(args.seqlengths);
	free(args.seqlist);
	free(fnames);
	running--;
	return NULL;
    }

    args.root = (ll_item *) malloc(sizeof(ll_item));
    if (args.root == NULL)
    {
	fastq_close(args.fastq);
	free(args.seqlist);
	free(args.seqlengths);
	free(fnames);
	running--;
	return PyErr_NoMemory();
    }
    args.root->fpos = -1L;
    args.root->next = NULL;

    args.pyhitseqs = PyList_New(0);

    init_stats(args.seqlist);

    profiles_n = 0;

    // start threads {{{3

    threads = (pthread_t *) malloc(sizeof(pthread_t) *nthreads);

    Py_BEGIN_ALLOW_THREADS

    for(threadi=0; threadi<nthreads; threadi++)
    {

	err = pthread_create(threads+threadi, NULL, 
		(void *(*)(void *)) scan_filepart,
		(void *) &args);

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

	// DBG("created thread #%d (%li)", threadi, pthread2long(threads[threadi]));
    }

    ret = NULL; // will indicate error if not set

    // end threads, get results {{{3

    if (threadi != -1)
	for(threadi=0; threadi<nthreads; threadi++)
	    pthread_join(threads[threadi], NULL);

    Py_END_ALLOW_THREADS

    if (threadi != -1)
    {

	if (exception == NULL)
	{
	    // convert return value

	    for(i=0,item=args.root->next; item!=NULL; i++)
		item = item->next;

	    hits = PyTuple_New(i);
	    for(i=0,item=args.root->next; item!=NULL; i++,item=item->next)
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

	    pystats = engine_stats(self, Py_BuildValue("()"));
	    ret = Py_BuildValue("{sOsOsO}",
		    "hits", hits,
		    "stats", pystats,
		    "hitseqs", args.pyhitseqs);
	    Py_XDECREF(hits);
	    Py_XDECREF(pystats);
	    Py_XDECREF(args.pyhitseqs);
	}
    }

    // clean up {{{3

    free(threads);
    free_ll(args.root);
    fastq_close(args.fastq);
    free(args.seqlist);
    free(args.seqlengths);

    if (exception != NULL) {
	PyErr_SetString(exception, errstr);
	exception = NULL;
    }

    profile_dump();

    running--;
    return ret;
}

/* engine.stop {{{2 */

    static PyObject *
engine_stop(PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, ""))
	return NULL;

    stop++;

    Py_RETURN_NONE;
}

/* engine.get_config {{{2 */

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

/* engine.config {{{2 */

    static PyObject *
engine_config(PyObject *self, PyObject *args, PyObject *kw)
{
    static char *kwl[] = { "maxerrors", "minoverlap", "minreadlength", "nthreads", "Amin", "Azero", NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kw, 
		"|iiiicc", kwl, &maxerrors, &minoverlap, &minreadlength, &nthreads, &Amin, &Azero))
	return NULL;

    Py_RETURN_NONE;
}

/* engine.test {{{2 */

    static PyObject *
engine_test(PyObject *self, PyObject *args)
{
    profile_dump();

    Py_RETURN_NONE;
}


/* initialization {{{1 */

/* PyMethodDef {{{2 */

static PyMethodDef methods[] = {
    {"test",  engine_test, METH_VARARGS,
	"test() -- perform some test.\n"},
    {"config",  (PyCFunction)engine_config, METH_VARARGS | METH_KEYWORDS,
	"config(**kwargs) -- configure the engine.\n"
	"arguments:\n"
	"'maxerrors : maximum number of base mismatches in sequence alignment\n"
	"'minoverlap : minimum number of base overlap for beginning/end hits\n"
	"'minreadlength : ignore reads shorter than this\n"
	"'nthreads : number of threads to use for scanning\n"
	"'Amin : nucleotides with quality ASCII value lower than this are discarded\n"
	"'Azero : ASCII value that corresponds to Q=0 (depends on FastQ format)\n"},
    {"get_config", engine_get_config, METH_VARARGS,
	"get_config() -- get the current config as dictionary.\n"},
    {"findseqs", engine_findseqs, METH_VARARGS,
	"findseqs(fname, sequences) -- finds occurences of base sequences in fastq files.\n"
	"arguments:\n"
	"'fname' : filename of fastq file or sequence of filenames of fastq files\n"
	"'sequences' : list of sequences to look for\n\n"
	"returns a dictionary with:\n"
	"'hits' : tuple of kvarq.engine.Hit\n"
	"'stats' : is the same dict as returned by a call to stats()\n"
	"'hitseqs' : tuple of base sequences corresponding to 'hits'\n"},
    {"stop",  engine_stop, METH_VARARGS,
	"stop() -- stops the scanning process.\n"},
    {"stats", engine_stats, METH_VARARGS,
	"stats() -- get statistics during scanning process.\n"
	"returns a dict containing:\n"
	"'readlengths' : tuple of number of occurences when accessed by read length\n"
	"'progress' : current progress of findseqs() ranging 0..1\n"
	"'sigints' : how many <CTRL-C> were caught since beginning of scan\n"
	"'nseqbasehits' : sum(hit_length), indexed by sequence as given to findseqs()\n"
	"'nseqhits' : number of hits, indexed by sequence as given to findseqs()\n"
	"'records_parsed' : total number of records parsed\n"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


/* init module {{{2 */

// if called initX : "module does not defined init function (initc)"
// if called initd and setup.py calls extensions 'd' : "dynamic module not initialized properly"
    PyMODINIT_FUNC
initengine(void)
{
    PyObject *mod, *obj;

    (void) Py_InitModule("engine", methods);

    // create object kvarq.engine.Hit

    mod = PyImport_ImportModule("collections");
    obj = PyObject_GetAttrString(mod, "namedtuple");
    Py_DECREF(mod);

    hittuple = PyObject_CallObject(obj,
	    Py_BuildValue("(sO)", "Hit", 
		Py_BuildValue("(sssss)",
		    "seq_nr",
		    "file_pos",
		    "seq_pos",
		    "length",
		    "readlength")));
    PyObject_SetAttrString(hittuple,"__doc__",PyString_FromString(
		"seq_nr : refers to the list of sequences in call to engine.findseqs\n"
		"file_pos : beginning of read (within decompressed data)\n"
		"seq_pos : places the beginning of the read relative to the beginning of the sequence (<0 if read overlaps only with beginning of sequence or read contains whole sequence; >0 if read overlaps only with end of sequence or read is contained within sequence)\n"
		"length : gives the number of overlapping basepairs\n"
		"readlength : length of the (quality trimmed) read containing the hit\n"
		));
    Py_DECREF(obj);
    engine_mod = PyImport_ImportModule("kvarq.engine");
    PyObject_SetAttrString(engine_mod, "Hit", hittuple);

    mod = PyImport_ImportModule("kvarq.fastq");
    fastq_exception = PyObject_GetAttrString(mod, "FastqFileFormatException");
    Py_DECREF(mod);
    exception = NULL;

    // initialize logging

    mod = PyImport_ImportModule("kvarq.log");
    obj = PyObject_GetAttrString(mod, "lo");
    lo_log = PyObject_GetAttrString(obj, "log");
    Py_DECREF(obj);
    Py_DECREF(mod);

    mod = PyImport_ImportModule("logging");
    obj = PyObject_GetAttrString(mod, "DEBUG");
    LOG_DEBUG = PyInt_AS_LONG(obj);
    Py_DECREF(obj);
    obj = PyObject_GetAttrString(mod, "INFO");
    LOG_INFO = PyInt_AS_LONG(obj);
    Py_DECREF(obj);
    obj = PyObject_GetAttrString(mod, "WARNING");
    LOG_WARNING = PyInt_AS_LONG(obj);
    Py_DECREF(obj);
    obj = PyObject_GetAttrString(mod, "ERROR");
    LOG_ERROR = PyInt_AS_LONG(obj);
    Py_DECREF(obj);
    obj = PyObject_GetAttrString(mod, "FATAL");
    LOG_FATAL = PyInt_AS_LONG(obj);
    Py_DECREF(obj);
    Py_DECREF(mod);

#ifdef _WIN32 // 32/64 bit windows
    SetConsoleCtrlHandler((PHANDLER_ROUTINE) CtrlHandler, TRUE);
#else // POSIX
    signal(SIGINT, sigint_cb);
#endif
}

// vim:ts=8:noet:sw=4:

