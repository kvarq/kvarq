import math

from kvarq.log import lo, tictoc

'''
from http://en.wikipedia.org/wiki/FASTQ_format#Encoding:

  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  ..LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                
                           -5....0........9.............................40 
                                 0........9.............................40 
                                    3.....9.............................40 
  0.2......................26...31........41                               

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
     (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
'''


class FastqFileFormatException(Exception):
    pass

class Fastq:


    ASCII = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ' + \
            '[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

    vendor_variants = (
            ('Sanger', range(0, 50), 0),
            ('Solexa', range(-5, 41), 31),
            ('Illumina 1.3+', range(0, 41), 31),
            ('Illumina 1.5+', range(3, 42), 31),
            ('Illumina 1.8+', range(0, 42), 0)
        )

    def __init__(self, fname, variant = None, dQ = None, fd=None):
        self.fname = fname

        if fd:
            self.fd = fd
        else:
            self.fd = open(self.fname, 'rb')
        self.fd.seek(0,2)
        self.size = self.fd.tell()

        min_pos, max_pos = self.min_max_score_check_file()
        lo.debug('min_pos=%d max_pos=%d' % (min_pos, max_pos))

        self.variants = []
        self.dQ = None
        for vendor_variant in self.vendor_variants:

            name, pos_range, dQ = vendor_variant

            if (min_pos - dQ) in pos_range and (max_pos - dQ) in pos_range:
                if self.dQ is None:
                    self.dQ = dQ
                self.variants.append(name)
                if self.dQ != dQ:
                    raise FastqFileFormatException(
                            'PHRED scores are ambiguous : ' + str(self.variants))
                lo.debug('compatible with ' + name)

        if self.dQ is None:
            raise FastqFileFormatException(
                    'PHRED score (%d..%d) incompatible with known formats' %
                    (min_pos, max_pos))

        self.Azero = self.ASCII[self.dQ]

        self.fd.seek(0)
        lines = [self.fd.readline() for i in range(4)]
        self.readlength = len(lines[1].strip('\r\n'))
        self.records_approx = self.size / sum([len(line) for line in lines])

        lo.info('readlength=%d records_approx=%d dQ=%d variants=%s' % (
                self.readlength, self.records_approx, self.dQ, str(self.variants)))


    def min_max_score_check_file(self, n=1000, points=10):
        ret_min = +999
        ret_max = -999
        for point in range(points):
            # (oversamples small files)
            self.fd.seek(self.size*point/points)
            if point > 0:
                self.seekback()
            while n > (points - 1 - point)*n/points:
                identifier = self.fd.readline().rstrip('\n\r')
                if not identifier: break
                if not identifier[0] == '@':
                    raise FastqFileFormatException(
                        'identifier (1st line of record) must begin with "@"')
                bases = self.fd.readline().rstrip('\n\r')
                if not set(bases).issubset(set('AGCTN')):
                    raise FastqFileFormatException(
                        'bases (2nd line of record) must contain only AGCTN')
                plus = self.fd.readline().rstrip('\n\r')
                if not (plus == '+' or (plus[0]=='+' and plus[1:] == identifier[1:])):
                    raise FastqFileFormatException(
                        'separator (3rd line of record) must be == "+" or "+(ident)"')
                phredstr = self.fd.readline().rstrip('\n\r')
                if not (len(bases) == len(phredstr) or (
                        len(bases) == len(phredstr)-1 and phredstr[-1] == '!' )):
                    raise FastqFileFormatException(
                        'bases must be ~ same length as phred score (2nd, 4th line)')
                try:
                    ret_min = min(ret_min, *[self.ASCII.index(x) for x in phredstr])
                    ret_max = max(ret_max, *[self.ASCII.index(x) for x in phredstr])
                except ValueError, e:
                    raise FastqFileFormatException(
                        'phred score (4th line of record) must contain only "%s"'%
                        self.ASCII)
                n -= 1

            if not identifier: break

        if not identifier:
            while self.fd.tell()<self.size:
                empty = self.fd.readline().rstrip('\n\r')
                if not empty == '':
                    raise FastqFileFormatException('non-empty line after empty line (fpos=%d'%
                            self.fd.tell())

        return ret_min, ret_max


    def A2Q(self, A):
        ''' translate PHRED ASCII value to Q value '''
        return self.ASCII.index(A) - self.dQ

    def Q2A(self, Q):
        ''' inverse of A2Q() '''
        return self.ASCII[Q + self.dQ]


    def Q2p(self, Q):
        ''' translate PHRED Q value to probability '''
        p = 10**(-.1*Q)
        # solexa prior to v1.3 should : p /= p+1
        return p


    def p2Q(self, p):
        ''' inverse of Q2p() '''
        # solexa prior to v1.3 should ...
        return int(-10 * math.log(p)/math.log(10))


    def Qs(self, n=1000):
        fd = open(self.fname)
        self.fd.seek(0)
        Qs = []
        for i in range(n):
            ident, seq, plus, scores = (self.fd.readline().strip()
                    for j in range(4))
            if not ident: break
            Qs.extend([self.A2Q(score) for score in scores])
        return Qs


    def lengths(self, Amin, n=1000, points=10):
        lengths = []
        for point in range(points):
            # (oversamples small files)
            self.fd.seek(self.size*point/points)
            if point > 0:
                self.seekback()

            while n > (points - 1 - point)*n/points:
                ident, seq, plus, scores = (self.fd.readline().strip()
                        for j in range(4))
                pos, length = self.cutoff(scores, Amin)
                if length>=0:
                    lengths.append(length)

                n -= 1
        return lengths

    def cutoff(self, scores, Amin):
        ''' returns ``pos, length`` of sequence that has >= ``Amin`` quality
            (use :py:meth:`Q2A` to convert quality value to ``ASCII``) '''
        length = -1
        pos_ = pos = 0
        for j, A in enumerate(scores):
            if ord(A) >= ord(Amin):
                if pos < 0: pos = j
            else:
                if pos >= 0 and length<j-pos:
                    length = j-pos
                    pos_ = pos
                pos = -1
        return pos_, length


    def readhit(self, hit):
        ''' :param hit: a :py:class:`kvarq.engine.Hit`
            :returns: a string base sequence '''
        if hit.seq_pos < 0:
            self.fd.seek(hit.file_pos-hit.seq_pos)
            return self.fd.read(hit.length)
        else:
            self.fd.seek(hit.file_pos)
            return self.fd.read(hit.length)

    def lineup(self):
        ''' sets file position to beginning of current line (or beginning of
            last line if it's already set to beginning of current line) '''
        pos = self.fd.tell()
        c = None
        while c!='\n' and pos>0:
            pos -= 1
            self.fd.seek(max(0, pos-1))
            c = self.fd.read(1)
        if pos == 0:
            self.fd.seek(0)

    def seekback(self):
        ''' moves file pointer to beginning of current record '''
        l = pos = None
        while pos!=0:
            # go one line up
            self.lineup()
            l = self.fd.readline(); self.lineup()
            if l[0] == '+':
                # go one line further up
                self.lineup()
                # previous + could be from quality score...
                l = self.fd.readline(); self.lineup()
                if l[0] == '+':
                    # ... if so, go one line further up
                    self.lineup()
                # go to beginning of record
                self.lineup()
                break
            pos = self.fd.tell()

    def dumpat(self, pos, Amin=None):
        ''' dumps the record at file position ``pos`` -- if ``Amin`` is specified
            then it also prints pos/length of sequence with given quality cutoff
            (use :py:meth:`Q2A` to convert quality value to ``ASCII``) '''
        self.fd.seek(pos)
        self.seekback()
        ident, seq, plus, scores =  self.readrecord()
        print ident
        print '-'.join(['%s(%d)'%(base, self.A2Q(A))
                for base, A in zip(seq, scores)])

        if Amin:
            pos, length = self.cutoff(scores, Amin)
            print 'Amin=%s -> pos=%d length=%d'%(Amin, pos, length)

    def readrecord(self):
        ''' reads one record; :py:attr:`fd` must point at first character of a
            record (use :py:meth:`seekback` first) '''
        ident, seq, plus, scores = (self.fd.readline().strip()
                for j in range(4))
        return ident, seq, plus, scores


    @tictoc('fastq.readhits')
    def readhits(self, hits):
        return [self.readhit(hit) for hit in hits]


if __name__ == '__main__':

    import argparse
    import json
    import logging

    from kvarq.util import TextHist

    parser = argparse.ArgumentParser(description='perform simple analysis on .fastq files')

    parser.add_argument('-v', '--verbose', action='count',
            help='show more generous information (on stdout)')
    parser.add_argument('-d', '--debug', action='store_true',
            help='output log information at a debug level')

    parser.add_argument('-n', '--reads', type=int, default=10000,
            help='how many reads to use for analysis')
    parser.add_argument('-m', '--parts', type=int, default=10,
            help='set to "1" to analyze all reads in beginning of file; set to "n" to distribute reads evenly across the whole file (more representative, but slower)')

    parser.add_argument('-Q', '--quality', type=int,
            help='show histogram of readlengths with specified quality score')

    parser.add_argument('fastq', help='name of .fastq file to analyze')
    
    args = parser.parse_args()
    fastq = args.fastq
    if fastq.endswith('.json'):
        lo.info('parsing ' + fastq)
        data = json.load(file(fastq))
        fastq = data['info']['fastq']
        lo.info('analyzing ' + fastq)
    fastq = Fastq(args.fastq)
    if args.debug:
        lo.setLevel(logging.DEBUG)

    if args.quality:
        Amin = fastq.Q2A(args.quality)
        n = args.reads
        m = args.parts
        lengths = []
        lo.debug('Q=%d p=%.5f Amin=%s'%(
                args.quality, fastq.Q2p(args.quality), Amin))
        for i in range(m):
            lengths += fastq.lengths(Amin, n/m, fastq.size*i/m)
        print TextHist().draw(sorted(lengths))

