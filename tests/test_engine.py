from kvarq import engine
from kvarq.fastq import Fastq, FastqFileFormatException
from kvarq import genes

import unittest
import os.path
import math
import random
import gzip
import tempfile


class FastqGenerator:

    ''' generate .fastq files that contain a specified set of sequences
        at a given quality score and a lot of random data '''

    def __init__(self, fname, force=False, variant='Sanger', debug=False):
        '''
        :param fname: name of output file
        :param force: whether to overwrite output file
        :param variant: a ``.fastq`` variant (see
            :py:attr:`kvarq.fastq.Fastq.vendor_variants)
        :param debug: whether to write debug information to stdout
        '''
        assert force or not os.path.exists(fname)
        self.fname = fname
        self.fd = file(fname, 'w')
        self.i = 0

        self.debug = debug

        self.variant = None
        for name, pos_range, dQ in Fastq.vendor_variants:
            if variant == name:
                self.variant = name
                self.dQ = dQ
                self.pos_min = pos_range[0]
                self.pos_max = pos_range[-1]
                break

        assert self.variant

    def write_record(self, bases, phredstr):
        ''' add record with given bases/PHRED score '''
        assert len(bases) == len(phredstr)
        self.fd.write('@GENERATED%010d\n'%self.i)
        self.i += 1
        self.fd.write(bases + '\n+\n' + phredstr+'\n')

    def size(self):
        ''' returns filesize '''
        return self.fd.tell()

    def p2A(self, p):
        ''' convert ``p`` value to ``ASCII`` value for this ``.fastq`` variant '''
        Q = int(-10 * math.log(p)/math.log(10))
        return Fastq.ASCII[Q + self.dQ]

    def randseq(self, n):
        ''' return a random base sequence of length ``n`` '''
        return ''.join(['ACGT'[random.randint(0,3)] for x in range(n)])

    def write_seq(self, seq, pmax=.05, left=0, right=0):
        ''' write a sequence and pad with random data on either side

            :param seq: base sequence to write
            :param pmax: maximum ``p`` value of bases in specified sequence
            :param left: number of random bases to add before ``seq``
                (with ``p>pmax``) 
            :param right: number of random bases to add after ``seq``
                (with ``p>pmax``) '''
        A = Fastq.ASCII
        Amin = self.p2A(pmax)
        Aidx = A.index(Amin)
        bases = self.randseq(left)
        bases+= seq
        bases+= self.randseq(right)
        phredstr = ''.join([A[random.randint(self.pos_min,Aidx-1)] for x in range(left)])
        phredstr+= ''.join([A[random.randint(Aidx,self.pos_max)] for x in range(len(seq))])
        phredstr+= ''.join([A[random.randint(self.pos_min,Aidx-1)] for x in range(right)])
        self.write_record(bases, phredstr)

    def cover_seq(self, seq, minoverlap, readlength, pmax=.05, left=10, right=10):
        ''' writes a record that cover a (long) sequence with specified parameters

            :param seq: base sequence
            :param minoverlap: minimum overlap between generated reads and ``seq``
            :param readlength: length of reads to generate (is also maximum overlap)
            :param pmax: maximum ``p`` value of bases in overlap
            :param left: number of random bases to add before overlap
                (with ``p>pmax``) 
            :param right: number of random bases to add after overlap
                (with ``p>pmax``) '''
        assert readlength>=minoverlap and len(seq)>=minoverlap

        overlap = random.randint(minoverlap, readlength)

        if overlap>len(seq):
            l = random.randint(0, readlength-len(seq))
            r = readlength-len(seq)-l
            seqx = self.randseq(l) + seq + self.randseq(r)

        else:
            if random.random()<.5:
                seqx = self.randseq(readlength-overlap) + seq[:overlap]
            else:
                seqx = seq[-overlap:] + self.randseq(readlength-overlap)

        assert len(seqx) == readlength
        self.write_seq(seqx, pmax=pmax, left=left, right=right)

    def flush(self):
        self.fd.flush()


class TestEngine(unittest.TestCase):

    def tmpfname(self):
        return os.path.splitext(__file__)[0] + '_tmp.fastq'

    def setUp(self):
        self.fname = os.path.join(os.path.dirname(__file__), 'test_engine.fastq')
        self.fname_1 = os.path.join(os.path.dirname(__file__), 'test_engine_1.fastq')
        self.fname_2 = os.path.join(os.path.dirname(__file__), 'test_engine_2.fastq')
        engine.config(nthreads=1)

    def tearDown(self):
        if os.path.exists(self.tmpfname()):
            os.unlink(self.tmpfname())


    def test_findseqs(self, gz=False):
        ''' find specified sequences in handwritten .fastq file '''
        engine.config(maxerrors=0, minoverlap=1000, minreadlength=3, Amin='!')
        seqs = (
            "CCC", # "CCCC" should be counted 2x ...
            "TTTT",
            "TATATATA",
            "TGTAG", # at beginning
            "ATATT", # at end
            "GAGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT",
            "...NACTTCCTCTCTACTGGTGTCGGCGGTGAAAGAGCTTACGTACTCTTCGAT...",
        )
        fname = self.fname
        if gz:
            fname += '.gz'
        hits = engine.findseqs(fname, seqs)['hits']

        if gz:
            f = gzip.GzipFile(fname, 'rb')
        else:
            f = file(fname, 'rb')

        x = [0] * len(seqs)

        for hit in hits:
            x[hit.seq_nr] += 1

            seq = seqs[hit.seq_nr]
            if hit.seq_pos<0:
                f.seek(hit.file_pos-hit.seq_pos)
                bps = f.read(hit.length)
            else:
                f.seek(hit.file_pos)
                bps = f.read(hit.length)
                seq = seq[hit.seq_pos:hit.seq_pos+hit[3]]

            assert bps == seq

        assert x == [19,1,0,1,1,1,1]

    def test_gz(self):
        ''' test engine gzip functionality '''
        self.test_findseqs(gz=True)

    def test_paired(self):
        engine.config(maxerrors=0, minoverlap=1000, minreadlength=3, Amin='!')
        seqs = (
            "CCC", # "CCCC" should be counted 2x ...
            "TTTT",
            "TATATATA",
            "TGTAG", # at beginning
            "ATATT", # at end
            "GAGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT",
            "...NACTTCCTCTCTACTGGTGTCGGCGGTGAAAGAGCTTACGTACTCTTCGAT...",
        )

        ret = engine.findseqs(self.fname, seqs)
        ret_12 = engine.findseqs((self.fname_1, self.fname_2), seqs)

        assert ret == ret_12

    def test_maxerror(self):
        ''' test different values for ``maxerror`` config parameter '''
        engine.config(minreadlength=25, minoverlap=25, Amin='!')
        seqs = (
            #GAGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT
            "CAGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT",
            #^ : 1 error
            "CTGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT",
            #^^: 2 errors
            "CTCCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT",
            #^^^: 3 errors
        )

        for maxerrors in range(4):
            engine.config( maxerrors=maxerrors )
            hits = engine.findseqs(self.fname, seqs)['hits']
            assert len(hits) == maxerrors


    def test_minoverlap(self):
        ''' test different values for ``minoverlap`` config parameter '''
        seqs = (
            "TCGATGCGATCTGTCAAGTCGGTGGCGGTA...", # end of sequence + junk
            "TCGATGCGATCTG.CAAGTCGGTGGCGGTA...", # end of sequence + junk + 1 error
            "...NTGAACGTATCGCCTCGAGGGACTT", # junk + beginning of sequence
            "...NTGAACGTATCG.CTCGAGGGACTT", # junk + beginning of sequence + 1 error
        )

        engine.config(
                maxerrors=0,
                minreadlength=25,
                minoverlap=30,
                Amin='!'
            )
        ret = engine.findseqs(self.fname, seqs)
        hits = ret['hits']
        assert len(hits)==1 and hits[0].seq_nr==0 and hits[0].seq_pos<0

        engine.config(maxerrors=0, minoverlap=25)
        hits = engine.findseqs(self.fname, seqs)['hits']
        assert len(hits)==2
        for hit in hits:
            assert hit[0]!=3 or hit[2]>0

        engine.config(maxerrors=1, minoverlap=25)
        hits = engine.findseqs(self.fname, seqs)['hits']
        assert len(hits)==4


    def test_Amin(self):
        seqs = (
                "GGAG",
                "CCGAC",
            )
        engine.config(Amin='H', minreadlength=4)
        ret = engine.findseqs(self.fname, seqs)

        assert len(ret['hits']) == 1
        assert ret['stats']['readlengths'][5] == 3
        assert ret['stats']['readlengths'][4] == 5

        engine.config(Amin='G')
        ret = engine.findseqs(self.fname, seqs)
        assert len(ret['hits']) == 2


    def test_hits(self):
        fname = self.tmpfname()
        fq = FastqGenerator(self.tmpfname(), force=True)
        seq = fq.randseq(51)

        minoverlap = 25
        readlength = 100
        pmax = .05
        n = 100
        for i in range(n):
            fq.cover_seq(seq,
                    minoverlap=minoverlap,
                    readlength=readlength,
                    pmax=pmax)
        fq.flush()
        #print "\033[94mfilesize=%.2f MB\033[m" % (fq.size() / 1024. / 1024.)

        fq = Fastq(fname)

        engine.config(
                nthreads=3,
                Amin=fq.Q2A(fq.p2Q(pmax)),
                maxerrors=0,
                minreadlength=random.randint(minoverlap, readlength),
                minoverlap=minoverlap
            )
        ret = engine.findseqs(fq.fname, [seq])

        assert ret['stats']['readlengths'][readlength] == n
        assert len(ret['hits']) == n

        if 0:
            print('hits=%d'%len(ret['hits']))
            print('readlenghts='+', '.join(['%dx %dbp'%(n, idx)
                    for idx,n in enumerate(ret['stats']['readlengths']) if n]))

        seqx = ''.join([i%minoverlap!=0 and b or {'A':'C','C':'G','G':'T','T':'A'}[b]
                    for i,b in enumerate(seq)])
        ret = engine.findseqs(fq.fname, [seqx])

        if 0:
            print('0123456789'*6)
            print(('*'+' '*(minoverlap-1))*6)
            print(seq)
            print(seqx)
            print(str(ret['hits']))

        assert ret['stats']['readlengths'][readlength] == n
        assert len(ret['hits']) == 0


    def test_fastq(self):

        file(self.tmpfname(), 'w').write('''_IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
#############################################
''')
        try:
            engine.findseqs(self.tmpfname(), [])
            assert False, "malformed @IDENTIFIER must raise FastqFileFormatException"
        except FastqFileFormatException:
            pass

        file(self.tmpfname(), 'w').write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
-
#############################################
''')
        try:
            engine.findseqs(self.tmpfname(), [])
            assert False, "malformed 3rd line must raise FastqFileFormatException"
        except FastqFileFormatException:
            pass


    def test_forward_fastq(self):
        engine.config(Amin='#', nthreads=2)
        for n in [3, 5, 7, 133]:
            for plus in ['+', '+IDENTIFIER']:
                for cr in ['\n', '\r\n']:
                    record = '@IDENTIFIER' + cr + 'A' * 80 + cr + \
                            plus + cr + '#' * 80 + cr
                    file(self.tmpfname(), 'wb').write(record * n)
                    Fastq(self.tmpfname())
                    ret = engine.findseqs(self.tmpfname(), ['A'*80])
                    assert len(ret['hits']) == n


if __name__ == '__main__': unittest.main()

