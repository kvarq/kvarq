
from kvarq.fastq import Fastq, FastqFileFormatException
from kvarq.log import lo

import unittest
import tempfile
import gzip
import os
import logging

from _util import lo_exceptor

class TestFastq(unittest.TestCase):

    def setUp(self):
        from kvarq.log import set_debug, set_warning
        #set_debug()
        #set_warning()
        self.gz = False
        self.tfastq = __file__ + '.fastq'

    def tearDown(self):
        for gz in ['', '.gz']:
            if os.path.exists(self.tfastq + gz):
                os.unlink(self.tfastq + gz)

    def ntf_write_fastq(self, content, variant=None):
        ntfn = self.tfastq
        if self.gz:
            ntfn += '.gz'
            gzf = gzip.GzipFile(ntfn, 'w')
            gzf.write(content)
            gzf.close()
        else:
            ntf = file(ntfn, 'w')
            ntf.write(content)
            ntf.close()
        return Fastq(ntfn, variant=variant)

    def ntf_write_quality(self, quality, variant=None):
        return self.ntf_write_fastq(
                '@IDENTIFIER\n' +
                'A' * len(quality) + '\n' +
                '+\n' +
                quality + '\n', variant=variant)

    def test_fastq_variant(self):
        lo.setLevel(logging.WARNING)

        fq = self.ntf_write_quality('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ');
        assert fq.dQ == 0 and set(fq.variants) == set(['Illumina 1.8+', 'Sanger']), \
                'FastQ variant should be Illumina 1.8+ / Sanger (dQ=0)'

        fq = self.ntf_write_quality(';<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh')
        assert fq.dQ == 31 and fq.variants == ['Solexa'], \
                'FastQ variant should be Solexa (dQ=31)'

        fq = self.ntf_write_quality('OPQRSTUVWXYZ[\]^_`abcdefgh')
        assert fq.dQ == 31 and fq.variants == ['Solexa', 'Illumina 1.3+', 'Illumina 1.5+'], \
                "FastQ variant should be ['Solexa', 'Illumina 1.3+', 'Illumina 1.5+'] (dQ=31)"
        try:
            fq = self.ntf_write_quality(';<=>?@ABCDEFGHI;<=>?@ABCDEFGHI')
            raise Exception(
                    'Quality score ";<=>?@ABCDEFGHI" ambiguous, should not be decoded!')
        except FastqFileFormatException:
            pass

        # specify valid vendor variants
        fq = self.ntf_write_quality(';<=>?@ABCDEFGHI;<=>?@ABCDEFGHI',
                variant='Sanger')
        fq = self.ntf_write_quality(';<=>?@ABCDEFGHI;<=>?@ABCDEFGHI',
                variant='Solexa')

        # specify invalid vendor variant, must emit warning
        lo_assert = lo_exceptor('seems not to be compatible', logging.WARNING, True)
        fq = self.ntf_write_quality(';<=>?@ABCDEFGHI;<=>?@ABCDEFGHI',
                variant='Illumina 1.3+')
        lo_assert()

        lo.setLevel(logging.INFO)


    def test_fastq_format(self):

        try:
            fq = self.ntf_write_fastq('''IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
############################################
''')
            assert False, 'identifier without "@" must raise exception'
        except FastqFileFormatException:
            pass

        try:
            fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTAXGTACGTACGTACGTACGTACGT
+
############################################
''')
            assert False, 'bases must not contain "X"'
        except FastqFileFormatException:
            pass

        try:
            fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
#############################################
''')
            assert False, 'phred string must be same length as bases'
        except FastqFileFormatException:
            pass

        try:
            fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+text
############################################
''')
            assert False, 'bases/phred score must be separated by "+" only'
        except FastqFileFormatException:
            pass

        try:
            fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
############################################

@IDENTIFIER
''')
            assert False, 'there must be no non-empy lines after empty lines'
        except FastqFileFormatException:
            pass

    def test_gz(self):
        ''' repeats tests with gzipped fastq files '''
        self.gz = True
        self.test_fastq_format()
        self.test_fastq_variant()
        self.gz = False

if __name__ == '__main__': unittest.main()

