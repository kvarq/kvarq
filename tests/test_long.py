
from kvarq import genes
from kvarq import engine
from kvarq import analyse
from kvarq.util import ProgressBar

import unittest
import urllib2
import shutil
import gzip
import os, os.path
import sys
from kvarq.fastq import Fastq
import json


MTB98url = 'ftp://ftp.broad.mit.edu/pub/annotation/mtuberculosis/diversity/MTB_98_1833.fastq.gz'

class TestLong(unittest.TestCase):

    def setUp(self):

        engine.config(nthreads=8)
        self.fname = os.path.join(os.path.dirname(__file__), 'MTB_98_1833.fastq')

        if not os.path.isfile(self.fname):
            basename = os.path.basename(self.fname)
            print >>sys.stderr, 'cannot find {} -- attempt to download from web'.format(basename)

            r = urllib2.urlopen(MTB98url)

            gzfname = self.fname + '.gz'
            with open(gzfname, 'w') as gzf:
                print >>sys.stderr, 'downloading...'

                m = r.info()
                total = left = int(m.dict['content-length'])
                pb = ProgressBar(total)
                pb.start()
                while left>0:
                    buf = r.read(min(10240, left))
                    gzf.write(buf)
                    left -= len(buf)

                    pb.update(total - left)
                    print >>sys.stderr, str(pb)

                with open(self.fname, 'w') as f:
                    print >>sys.stderr, 'unpacking...'
                    gzf = gzip.GzipFile(gzfname)
                    shutil.copyfileobj(gzf, f)

            os.unlink(gzfname)


    def test_long(self):

        tests = genes.phylo.lineages_SNPs
        analyser = analyse.Analyser(Fastq(self.fname), tests)
        analyser.scan()
        analyser.update()

        print json.dumps(analyser.serialize())

