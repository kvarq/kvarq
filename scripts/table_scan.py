
import argparse
import os
import os.path
import sys
import logging
import subprocess
import time
import shlex
#import postmortem

try:
    from kvarq import cli
except ImportError:
    sys.stderr.write('could not find module kvarq -- make sure it is in the PYTHONPATH\n')
    sys.exit(-1)

from table_utils import csv_xls_wrapper

parser = argparse.ArgumentParser(
        description='runs kvarq CLI over a list of .fastq '
        'files specified in a .csv/.xls table and saves the resulting .json '
        'files as well as a kvarq.log file into an output directory and a '
        'table_scan.log in the local directory -- by default, testsuites will '
        'be taken from the directory "testsuites/" in the local directory (see '
        '--flags argument below to change this default behavior)')

defaultflags = '-l kvarq.log -t testsuites/ scan -p'
parser.add_argument('-f', '--flags', default=defaultflags,
        help='flags to be passed on to kvarq/cli.py (default="%s")'% defaultflags)

parser.add_argument('-i', '--directory', default='.',
        help='where to look for .fastq files (defaults to ".")')

logfn = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.log'
parser.add_argument('-l', '--log', default=logfn,
        help='where to save log (defaults to %s)'%logfn)
parser.add_argument('-d', '--debug', action='store_true',
        help='show debug information on console')
parser.add_argument('-o', '--save-output', action='store_true',
        help='save output of kvarq to .txt files')
parser.add_argument('-M', '--no-mtime', action='store_true',
        help='ignore modification time (normally, .fastq is not scanned again if newer .json file already exists in directory')

parser.add_argument('-H', '--no-header', action='store_true',
        help='do not discard first row of table')
parser.add_argument('-c', '--column', type=int, default=0,
        help='specify column containing fastq name (defaults to 0=first column)')

parser.add_argument('table',
        help='.csv/.xls input file containing list of .fastq files to scan')
parser.add_argument('output_directory',
        help='where to store .json files')

args = parser.parse_args()
kvarq_args = shlex.split(args.flags.format(logfn=logfn))


# set up logging
lo = logging.getLogger('table_scan')
ch = logging.StreamHandler()
fh = logging.FileHandler(args.log)
ft = logging.Formatter('[%(asctime)s] -%(levelname)s- l_%(lineno)d :: %(message)s')
ch.setFormatter(ft)
fh.setFormatter(ft)
lo.setLevel(logging.DEBUG)
if args.debug:
    ch.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.INFO)
fh.setLevel(logging.INFO)
lo.addHandler(ch)
lo.addHandler(fh)

lo.info('started : ' + str(sys.argv))


if not os.path.isdir(args.output_directory):
    os.mkdir(args.output_directory)

# src: http://stackoverflow.com/questions/6796492/python-temporarily-redirect-stdout-stderr
class RedirectStdStreams(object):
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

header = None
table = csv_xls_wrapper(args.table)
logfn = os.path.join(args.output_directory, 'kvarq.log')

# iterate over table
for row in table:

    if not header and not args.no_header:
        header = row
        continue

    # find .json
    fname = row[args.column]
    if os.path.exists(fname):
        # can be specified as path
        fastq = fname
        fname = os.path.basename(fname)
    else:
        # or only filename -> find recursively
        fname = os.path.basename(fname)
        if not fname.endswith('.fastq'):
            fname += '.fastq'
        lo.debug('looking for %s...'% fname)
        for dirpath, dirnames, filenames in os.walk(args.directory):
            if fname in filenames:
                break
        if not fname in filenames:
            lo.error('could not locate file : ' + fname)
            continue
        fastq = os.path.join(dirpath, fname)

    jsonfn = os.path.join(args.output_directory, os.path.splitext(fname)[0] + '.json')
    txt = os.path.join(args.output_directory, os.path.splitext(fname)[0] + '.txt')

    if not args.no_mtime and os.path.exists(jsonfn) and \
            os.path.getmtime(fastq)<os.path.getmtime(jsonfn):
        lo.info('not scanning %s because %s exists and is newer'%(fastq, jsonfn))
        continue

    # run kvarq.cli
    lo.info('scanning %s : kvarq.cli(%s)'%(fastq, kvarq_args + [fastq, jsonfn]))

    fstdout = None
    if args.save_output:
        fstdout = file(txt, 'w')

    try:
        with RedirectStdStreams(fstdout):
            cli.main(kvarq_args + [fastq, jsonfn])
    except Exception, e:
        lo.error('error while scanning %s : %s'%(jsonfn, str(e)))

