
import math
import gzip
import os.path
import collections

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

    VendorProperties = collections.namedtuple(
            'VendorProperties', ['Qrange', 'dQ'])

    vendor_variants = dict((
            ('Sanger', VendorProperties(range(0, 50), 0)),
            ('Solexa', VendorProperties(range(-5, 41), 31)),
            ('Illumina 1.3+', VendorProperties(range(0, 41), 31)),
            ('Illumina 1.5+', VendorProperties(range(3, 42), 31)),
            ('Illumina 1.8+', VendorProperties(range(0, 62), 0))
        ))

    def __init__(self, fname, variant=None, fd=None, paired=False, quiet=False):
        '''
        open ``.fastq`` or ``.fastq.gz`` file and determine its
        variant (setting attribute ``.Azero`` accordingly)

        :param fname: name of file to open
        :param variant: specify one of ``.vendor_variants`` -- if none
            is specified, then the PHRED score of the fastq file is
            analyzed and
        :param fd: specify a file descriptor to use instead of
            opening ``fname``
        :param paired: include second file in a paired set if it is
            available (i.e. specify "file_1.fastq" as input file and
            "file_2.fastq" will be included in functions ``.filesize()``
            and ``.filenames()``)
        '''
        self.fname = fname

        if fd:
            self.fd = fd
        else:
            self.fd = None

        if self.fname.endswith('.fastq.gz'):
            self.gz = True
            if not self.fd:
                self.fd = gzip.GzipFile(self.fname, 'rb')
        elif self.fname.endswith('.fastq'):
            self.gz = False
            if not self.fd:
                self.fd = open(self.fname, 'rb')
        else:
            raise FastqFileFormatException(
                        'fastq file must have extension ".fastq" or ".fastq.gz"')

        # save second name of base if exists
        self.fname2 = None
        if paired:
            base = fname[:fname.rindex('.fastq')]
            if base[-2:] == '_1':
                fname2 = base[:-2] + '_2' + fname[fname.rindex('.fastq'):]
                if os.path.exists(fname2):
                    lo.info('including paired file "%s"' % fname2)
                    self.fname2 = fname2

        if sum(self.filesizes()) == 0:
            raise FastqFileFormatException('cannot scan empty file')

        # scan some records
        min_pos, max_pos = self.min_max_score_check_file()
        lo.debug('min_pos=%d max_pos=%d' % (min_pos, max_pos))

        if variant and variant not in self.vendor_variants:
            raise FastqFileFormatException(
                    'unknown vendor variant "%s"' % variant)

        # create list of variants compatible with PHRED scores
        variants = []
        dQs = []
        for name, vendor_variant in Fastq.vendor_variants.items():

            if ((min_pos - vendor_variant.dQ) in vendor_variant.Qrange
                    and (max_pos - vendor_variant.dQ) in vendor_variant.Qrange):
                dQs.append(vendor_variant.dQ)
                variants.append(name)

        if variant is None:
            # set variant from guesses
            if not variants:
                raise FastqFileFormatException(
                        'could not find any suitable fastq vendor variant')
            if len(set(dQs)) > 1:
                raise FastqFileFormatException(
                        'cannot determine dQ with guessed vendor variants "%s"'
                        % str(variants))
            self.variants = variants
            self.dQ = dQs[0]
        else:
            # check specified variant
            if variant not in variants:
                lo.warning('specified vendor variant "%s" seems not to be '
                        'compatible with file' % variant)
            self.variants = [variant]
            self.dQ = self.vendor_variants[variant].dQ


        self.Azero = self.ASCII[self.dQ]

        # estimate readlength/records_approx
        self.fd.seek(0)
        lines = [self.fd.readline() for i in range(4)]
        self.readlength = len(lines[1].strip('\r\n'))
        if self.gz:
            self.records_approx = None
        else:
            self.records_approx = os.path.getsize(self.fname) / len(''.join(lines))
            if self.fname2 is not None:
                self.records_approx *= 2

        # output some infos
        if not quiet:
            if self.gz:
                lo.info('gzipped fastq : readlength=? records_approx=? dQ=%d variants=%s' % (
                        self.dQ, str(self.variants)))
            else:
                lo.info('fastq : readlength=%d records_approx=%d dQ=%d variants=%s' % (
                        self.readlength, self.records_approx, self.dQ, str(self.variants)))


    def filesizes(self):
        ''' returns list of filesize(s) -- see ``paired`` parameter in
            :py:meth:`.__init__` '''
        return [os.path.getsize(fname) for fname in self.filenames()]

    def filenames(self):
        ''' returns list of filename(s) -- see ``paired`` parameter in
            :py:meth:`.__init__` '''
        if self.fname2 is not None:
            return [self.fname, self.fname2]
        return [self.fname]

    def min_max_score_check_file(self, n=1000, points=10):
        '''
        check fastq file format and return min/max PHRED score values

        :param n: number of records to scan
        :param points: number of points within file to scan for records;
            this value is ignored for gzipped fastq files
        :returns: minimum and maximum value of PHRED score (index within
            ``ASCII``)
        '''
        ret_min = +999
        ret_max = -999
        self.fd.seek(0)

        if self.gz:
            lo.debug('gzipped fastq : scan %d points at start only' % n)

        for point in range(points):

            if not self.gz and point > 0:
                # (oversamples small files)
                self.fd.seek(os.path.getsize(self.fname)*point/points)
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
            while True:
                line = self.fd.readline()
                if not line: break
                if not line.rstrip('\r\n') == '':
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


    def lengths(self, Amin, n=1000, points=10):
        '''
        samples length of quality trimmed records

        :param Amin: minimum PHRED value
        :param n: number of records to sample
        :param points: number of points within file to scan for records;
            this value is ignored for gzipped fastq files
        :returns: list of quality trimmed record lengths ``n`` items
        '''
        self.fd.seek(0)

        if self.gz:
            lo.debug('gzipped fastq : scan %d points at start only' % n)

        lengths = []
        for point in range(points):

            if not self.gz and point > 0:
                self.fd.seek(os.path.getsize(self.fname)*point/points)
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
        ''' moves file pointer to beginning current (if after ``plus``) or previous
            (if before ``plus``) record '''
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
        print(ident)
        print('-'.join(['%s(%d)'%(base, self.A2Q(A))
                for base, A in zip(seq, scores)]))

        if Amin:
            pos, length = self.cutoff(scores, Amin)
            print('Amin=%s -> pos=%d length=%d'%(Amin, pos, length))

    def readrecord(self):
        ''' reads one record; :py:attr:`fd` must point at first character of a
            record (use :py:meth:`seekback` first) '''
        ident, seq, plus, scores = (self.fd.readline().strip()
                for j in range(4))
        return ident, seq, plus, scores

    def readrecordat(self, hit):
        ''' :param hit: a :py:class:`kvarq.engine.Hit`
            :returns: the four .fastq files representing the record '''
        self.fd.seek(hit.file_pos)
        self.seekback()
        ident, seq, plus, scores = self.readrecord() # previous record
        ident, seq, plus, scores = self.readrecord() # our record
        return '\n'.join([ident, seq, plus, scores]) + '\n'

    @tictoc('fastq.readhits')
    def readhits(self, hits):
        return [self.readhit(hit) for hit in hits]

