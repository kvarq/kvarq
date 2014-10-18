
import unittest, os.path, logging
from distutils.version import StrictVersion

from kvarq.testsuites import discover_testsuites, load_testsuites, update_testsuites
from kvarq.analyse import TestsuiteVersionConflictException
from kvarq.log import lo


here_dir = os.path.abspath(os.path.dirname(__file__))
testsuites_alt = os.path.join(here_dir, 'override_testsuites')


class TestsuiteTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        lo.setLevel(logging.WARNING)
        cls.testsuite_paths = discover_testsuites([testsuites_alt])
        cls.testsuites = load_testsuites(cls.testsuite_paths, ['MTBC'])

    @classmethod
    def tearDownClass(cls):
        lo.setLevel(logging.INFO)

    def test_update_testsuites(self):

        v = StrictVersion(self.testsuites['MTBC/test'].version)
        # load by full name
        testsuites = {}
        update_testsuites(testsuites,
                {'MTBC/test': str(v)},
                self.testsuite_paths
            )
        assert testsuites.keys() == ['MTBC/test']
        # load by short name
        update_testsuites(testsuites,
                {'test': str(v)},
                self.testsuite_paths
            )
        assert set(testsuites.keys()) == set(['MTBC/test', 'test'])
        assert testsuites['test'] == testsuites['MTBC/test']

        # load compatible
        vv = list(v.version)
        vv[1] -= 1
        v.version = vv
        update_testsuites(testsuites,
                {'test': str(v)},
                self.testsuite_paths
            )
        assert set(testsuites.keys()) == set(['MTBC/test', 'test'])

        # load incompatbile 1/2
        vv[1] += 2
        v.version = vv
        try:
            update_testsuites(testsuites,
                    {'test': str(v)},
                    self.testsuite_paths
                )
            assert False, 'future minor version specified; should fail'
        except TestsuiteVersionConflictException:
            pass

        # load incompatbile 1/2
        vv[1] -= 1
        vv[0] -= 1
        v.version = vv
        try:
            update_testsuites(testsuites,
                    {'test': str(v)},
                    self.testsuite_paths
                )
            assert False, 'different major version specified; should fail'
        except TestsuiteVersionConflictException:
            pass


if __name__ == '__main__': unittest.main()

