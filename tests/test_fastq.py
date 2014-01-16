
from kvarq.fastq import Fastq, FastqFileFormatException

import unittest
import tempfile


class TestFastq(unittest.TestCase):

    def setUp(self):
        from kvarq.log import set_debug, set_warning
        #set_debug()
        set_warning()

    def test_fastq_variant(self):
        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAA
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
''')
        ntf.flush()
        fq = Fastq(ntf.name, fd=ntf.file)
        assert fq.dQ == 0 and set(fq.variants) == set(['Illumina 1.8+', 'Sanger']), \
                'FastQ variant should be Illumina 1.8+ / Sanger (dQ=0)'

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACTCTA
+
;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
''')
        ntf.flush()
        fq = Fastq(ntf.name, fd=ntf.file)
        assert fq.dQ == 31 and fq.variants == ['Solexa'], \
                'FastQ variant should be Solexa (dQ=31)'

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTAC
+
OPQRSTUVWXYZ[\]^_`abcdefgh
''')
        ntf.flush()
        fq = Fastq(ntf.name, fd=ntf.file)
        assert fq.dQ == 31 and fq.variants == ['Solexa', 'Illumina 1.3+', 'Illumina 1.5+'], \
                "FastQ variant should be ['Solexa', 'Illumina 1.3+', 'Illumina 1.5+'] (dQ=31)"

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACCAGT
+
;<=>?@ABCDEFGHI;<=>?@ABCDEFGHI
''')
        ntf.flush()
        try:
            fq = Fastq(ntf.name, fd=ntf.file)
            raise Exception(
                    'Quality score ";<=>?@ABCDEFGHI" ambiguous, should not be decoded!')
        except FastqFileFormatException:
            pass


    def test_fastq_format(self):

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
############################################
''')
        ntf.flush()
        try:
            fq = Fastq(ntf.name, fd=ntf.file)
            assert False, 'identifier without "@" must raise exception'
        except FastqFileFormatException:
            pass

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTAXGTACGTACGTACGTACGTACGT
+
############################################
''')
        ntf.flush()
        try:
            fq = Fastq(ntf.name, fd=ntf.file)
            assert False, 'bases must not contain "X"'
        except FastqFileFormatException:
            pass

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
#############################################
''')
        ntf.flush()
        try:
            fq = Fastq(ntf.name, fd=ntf.file)
            assert False, 'phred string must be same length as bases'
        except FastqFileFormatException:
            pass

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+text
############################################
''')
        ntf.flush()
        try:
            fq = Fastq(ntf.name, fd=ntf.file)
            assert False, 'bases/phred score must be separated by "+" only'
        except FastqFileFormatException:
            pass

        ntf = tempfile.NamedTemporaryFile()
        ntf.write('''@IDENTIFIER
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
############################################

@IDENTIFIER
''')
        ntf.flush()
        try:
            fq = Fastq(ntf.name, fd=ntf.file)
            assert False, 'there must be no non-empy lines after empty lines'
        except FastqFileFormatException:
            pass

if __name__ == '__main__': unittest.main()
