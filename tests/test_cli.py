
import unittest
import sys, os.path, os, logging, time, tempfile
from cStringIO import StringIO

from kvarq import VERSION
import kvarq.cli
from kvarq.log import lo

here_dir = os.path.abspath(os.path.dirname(__file__))
MTBC_fastq1 = os.path.join(here_dir, 'fastqs', 'L3_N1014_hits_5k.fastq')
MTBC_fastq2 = os.path.join(here_dir, 'fastqs', 'N0116_1_hits_1k.fastq')
root_dir = os.path.join(here_dir, os.pardir)
testsuites_alt = os.path.join(here_dir, 'override_testsuites')


class CliTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cwd = os.getcwd()
        os.chdir(here_dir)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.cwd)


    def main(self, args, err=0):

        stdout = sys.stdout
        stderr = sys.stderr
        strout = sys.stdout = StringIO()
        strerr = sys.stderr = StringIO()

        try:
            kvarq.cli.main(args)
        except SystemExit as e:
            sys.stdout = stdout
            sys.stderr = stderr
            assert e.code == err, (
                    'expected error code=%d' % err
                    + '\n\nstdout: ' + strout.getvalue()
                    + '\n\nstderr: ' + strerr.getvalue()
                )

        return strout.getvalue(), strerr.getvalue()


    def test_version(self):
        out, err = self.main(['version'])
        assert out.strip('\n\r') == VERSION


    def test_usage(self):

        out, err = self.main(['-h'])
        assert out.startswith('usage:'), 'expected usage output'


    def test_load_testsuites(self):


        def get_testsuites(out):
            line = [line for line in out.split('\n')
                    if line.startswith('testsuites=')][0]
            return set([testsuite
                for testsuite in line[line.index('=')+1:].split(',')
                if testsuite])

        def name_and_version(testsuites):
            return set([testsuite[:testsuite.find('[')]
                for testsuite in testsuites])

        def name_only(testsuites):
            return set([testsuite[:testsuite.find('-')]
                for testsuite in testsuites])

        cwd = os.getcwd()

        out, err = self.main(['info'])
        assert get_testsuites(out) == set()

        # suppress "loaded testsuite" messages
        lo.setLevel(logging.WARNING)

        MTBC_testsuites = set(['MTBC/phylo', 'MTBC/resistance', 'MTBC/spoligo'])

        # select testsuite
        out, err = self.main(['info', '-l', 'MTBC/phylo'])
        assert name_only(get_testsuites(out)) == set(['MTBC/phylo'])

        # select group
        out, err = self.main(['info', '-l', 'MTBC'])
        assert name_only(get_testsuites(out)) == MTBC_testsuites

        # select testsuite by filename
        path = os.path.join(testsuites_alt, 'MTBC', 'phylo.py')
        out, err = self.main(['info', '-l', path])
        assert name_and_version(get_testsuites(out)) == set(['MTBC/phylo-0.0'])

        # override testsuite directory using switch
        out, err = self.main(['info', '-l', 'MTBC/phylo'])
        assert name_and_version(get_testsuites(out)) != set(['MTBC/phylo-0.0'])
        out, err = self.main(['-t', testsuites_alt, 'info', '-l', 'MTBC/phylo'])
        assert name_and_version(get_testsuites(out)) == set(['MTBC/phylo-0.0'])

        # override testsuite directory using KVARQ_TESTSUITES
        os.environ['KVARQ_TESTSUITES'] = testsuites_alt
        out, err = self.main(['info', '-l', 'MTBC/phylo'])
        assert name_and_version(get_testsuites(out)) == set(['MTBC/phylo-0.0'])
        del os.environ['KVARQ_TESTSUITES']

        # time load all
        t0 = time.time()
        out, err = self.main(['info', '-L'])
        assert len(name_only(get_testsuites(out))) > 4
        if time.time() - t0 > 2:
            lo.warning('loading all testsuites takes %.2f' % (time.time() - t0))

        lo.setLevel(logging.INFO)
        os.chdir(cwd)


    def scan_illustrate(self, MTBC_fastq, scan_params=[]):

        ntf = tempfile.NamedTemporaryFile(delete=False)
        lo.setLevel(logging.WARNING)

        try:

            t0 = time.time()
            out, err = self.main(['scan', '-l', 'MTBC', '-f'] + scan_params + [
                MTBC_fastq, ntf.name])

            if time.time() - t0 > 10:
                lo.warning('scanning of %s took %.2fs' % (
                    os.path.basename(MTBC_fastq), time.time() - t0))

            out, err = self.main(['illustrate', '-r',
                ntf.name])

        finally:
            lo.setLevel(logging.INFO)
            ntf.close()
            os.remove(ntf.name)

        return out, err

    def test_scan_illustrate(self):

        lo.setLevel(logging.WARNING)

        # do a complete scan with a known result; specify variatn
        out, err = self.scan_illustrate(MTBC_fastq1,
                ['--variant', 'Illumina 1.8+'])
        for resistance in [
                'Streptomycin resistance::SNP781687AG=rpsL.K43R',
                'Ethambutol resistance::SNP4247431GT=embB.M306I',
                'Isoniazid resistance [2155168CG=katG.S315T]',
                'Rifampicin resistance (RRDR) [761139CG=rpoB.H445D 761140AG=rpoB.H445R]',
                'remark: low coverage (RRDR below 10x)']:
            assert resistance in out, MTBC_fastq1 + ' should have ' + resistance
        assert 'lineage 3' in out, MTBC_fastq1 + ' should be Lineage 3'

        # another tested .fsatq
        out, err = self.scan_illustrate(MTBC_fastq2)
        for resistance in [
                'Streptomycin resistance::SNP781687AG=rpsL.K43R',
                'remark: low coverage (RRDR below 10x)']:
            assert resistance in out, MTBC_fastq2 + ' should have ' + resistance
        assert 'lineage 2' in out, MTBC_fastq2 + ' should be Lineage 2'

        lo.setLevel(logging.INFO)


if __name__ == '__main__': unittest.main()

