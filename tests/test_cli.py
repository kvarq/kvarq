
import unittest
import sys

import kvarq.cli


class CliTest(unittest.TestCase):

    def test_real(self):

        sys.argv = ['','-h']
        print 'ARGV='+str(sys.argv)

        try:
            kvarq.cli.main()
        except SystemExit, e:
            pass

if __name__ == '__main__': unittest.main()
