from kvarq import engine
from kvarq.fastq import Fastq, FastqFileFormatException
from kvarq import genes

import unittest
import os.path
import math
import random


class FastqGenerator:

    def __init__(self, fname, force=False, variant='Sanger'):
        assert force or not os.path.exists(fname)
        self.fname = fname
        self.fd = file(fname, 'w')
        self.i = 0

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
        assert len(bases) == len(phredstr)
        self.fd.write('@GENERATED%010d\n'%self.i)
        self.i += 1
        self.fd.write(bases + '\n+\n' + phredstr+'\n')

        if 0:
            print bases
            print phredstr
            print ''.join([Fastq.ASCII.index(a)>=13+self.dQ and '*' or '_' for a in phredstr])

    def p2A(self, p):
        Q = int(-10 * math.log(p)/math.log(10))
        return Fastq.ASCII[Q + self.dQ]

    def randseq(self, n):
        return ''.join(['ACGT'[random.randint(0,3)] for x in range(n)])

    def write_seq(self, seq, pmax=.05, left=0, right=0):
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

    def cover_seq(self, seq, n, minoverlap, readlength, pmax=.05, left=10, right=10):
        assert readlength>=minoverlap and len(seq)>=minoverlap

        for i in range(n):
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

        self.flush()

    def flush(self):
        self.fd.flush()


class TestEngine(unittest.TestCase):

    def tmpfname(self):
        return os.path.splitext(__file__)[0] + '_tmp.fastq'

    def setUp(self):
        self.fname = os.path.join(os.path.dirname(__file__),'test_engine.fastq')
        engine.config(nthreads=1)
        if os.path.exists(self.tmpfname()):
            os.unlink(self.tmpfname())
        from kvarq.log import set_warning
        set_warning()


    def test_findseqs(self):
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
        hits = engine.findseqs(self.fname, seqs)['hits']

        f = file(self.fname, 'rb')

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


    def test_fuzzy(self):
        engine.config(minoverlap=1000, Amin='!')
        seqs = (
            "CAGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT", # 1 error
            "CTGCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT", # 2 errors
            "CTCCATGTGGAGCAACTTGTGGGAGCGCCGGGCAACGCCCTGTCTCTTAT", # 3 errors
        )

        for maxerrors in range(4):
            engine.config( maxerrors=maxerrors )
            hits = engine.findseqs(self.fname, seqs)['hits']
            assert len(hits) == maxerrors


    def test_overlap(self):
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
        fname = os.path.splitext(__file__)[0] + '_generated.fastq'
        fq = FastqGenerator(fname, force=True)
        seq = fq.randseq(51)

        minoverlap = 25
        readlength = 100
        pmax = .05
        n = 100
        fq.cover_seq(seq, n,
                minoverlap=minoverlap,
                readlength=readlength,
                pmax=pmax)

        fq = Fastq(fname)

        engine.config(
                nthreads=1,
                Amin=fq.Q2A(fq.p2Q(pmax)),
                maxerrors=0,
                minreadlength=random.randint(minoverlap, readlength),
                minoverlap=minoverlap
            )
        ret = engine.findseqs(fq.fname, [seq])

        assert ret['stats']['readlengths'][readlength] == n
        assert len(ret['hits']) == n

        if 0:
            print 'hits=%d'%len(ret['hits'])
            print 'readlenghts='+', '.join(['%dx %dbp'%(n, idx)
                    for idx,n in enumerate(ret['stats']['readlengths']) if n])

        seqx = ''.join([i%minoverlap!=0 and b or {'A':'C','C':'G','G':'T','T':'A'}[b]
                    for i,b in enumerate(seq)])
        ret = engine.findseqs(fq.fname, [seqx])

        if 0:
            print '0123456789'*6
            print ('*'+' '*(minoverlap-1))*6
            print seq
            print seqx
            print str(ret['hits'])

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
