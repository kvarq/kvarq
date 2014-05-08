
from kvarq.fastq import Fastq, FastqFileFormatException

import unittest
import tempfile
import gzip
import os


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

    def ntf_write_fastq(self, content):
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
        return Fastq(ntfn)

    def test_fastq_variant(self):
        fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAA
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
''')
        assert fq.dQ == 0 and set(fq.variants) == set(['Illumina 1.8+', 'Sanger']), \
                'FastQ variant should be Illumina 1.8+ / Sanger (dQ=0)'

        fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACTCTA
+
;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
''')
        assert fq.dQ == 31 and fq.variants == ['Solexa'], \
                'FastQ variant should be Solexa (dQ=31)'

        fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTAC
+
OPQRSTUVWXYZ[\]^_`abcdefgh
''')
        assert fq.dQ == 31 and fq.variants == ['Solexa', 'Illumina 1.3+', 'Illumina 1.5+'], \
                "FastQ variant should be ['Solexa', 'Illumina 1.3+', 'Illumina 1.5+'] (dQ=31)"

        try:
            fq = self.ntf_write_fastq('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACCAGT
+
;<=>?@ABCDEFGHI;<=>?@ABCDEFGHI
''')
            raise Exception(
                    'Quality score ";<=>?@ABCDEFGHI" ambiguous, should not be decoded!')
        except FastqFileFormatException:
            pass


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

