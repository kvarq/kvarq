'''
KvarQ command line client
'''

from kvarq import VERSION
from kvarq import genes
from kvarq import engine
from kvarq import analyse
from kvarq.util import ProgressBar, TextHist, json_dump
from kvarq.fastq import Fastq, FastqFileFormatException
from kvarq.log import lo, appendlog, set_debug, set_warning, format_traceback
from kvarq.config import default_config

import argparse
import sys
import threading
import time
import json
import os.path
import codecs
from pprint import pprint
import glob

ERROR_COMMAND_LINE_SWITCH = -1
ERROR_FASTQ_FORMAT_ERROR = -2
ERROR_JSON_EXISTS = -3

# utils {{{1

def traceit(type, value, tb):
    if hasattr(sys, 'ps1') or not sys.stderr.isatty():
        sys.__excepthook__(type, value, tb) # default hook
    else:
        import traceback, pdb
        traceback.print_exception(type, value, tb)
        print
        pdb.post_mortem(tb)

def heardEnter():
    # src : http://stackoverflow.com/questions/292095/polling-the-keyboard-in-python
    try:
        import select
        i,o,e = select.select([sys.stdin],[],[],0.0001)
        for s in i:
            if s == sys.stdin:
                input = sys.stdin.readline()
                return True
        return False
    except:
        # e.g. windows fails
        return False

def invalid_args(msg):
    print '\n*** {}\n\n'.format(msg)
    parser.print_help()
    sys.exit(ERROR_COMMAND_LINE_SWITCH)


# scan {{{1

def scan(args, testsuites):

    # prepare scanning {{{2

    try:
        fastq = Fastq(args.fastq)
    except FastqFileFormatException, e:
        lo.error('cannot open file %s : %s'%(args.fastq, str(e)))
        sys.exit(ERROR_FASTQ_FORMAT_ERROR)

    engine.config(
            nthreads=args.threads,
            maxerrors=args.errors,
            Amin=fastq.Q2A(args.quality),
            Azero=fastq.Azero,
            minreadlength=args.readlength,
            minoverlap=args.overlap
        )

    analyser = analyse.Analyser()

    if not args.force and os.path.exists(args.json):
        lo.error('will not overwrite file ' + args.json)
        sys.exit(ERROR_JSON_EXISTS)

    # do scanning {{{2

    mb = os.path.getsize(args.fastq) / 1024 / 1024
    lo.info('scanning {} ({} MB)...'.format(args.fastq, mb))
    t0 = time.time()

    class AnalyseThread(threading.Thread):

        def __init__(self, analyser):
            super(AnalyseThread, self).__init__(name='analyse-thread')
            self.analyser = analyser
            self.finished = False
            self.exception = None
            self.traceback = None

        def run(self):
            try:
                self.analyser.scan(fastq, testsuites, do_reverse=not args.no_reverse)
                self.finished = True
            except Exception, e:
                self.exception = e
                self.traceback = format_traceback(sys.exc_info())

    at = AnalyseThread(analyser)

    at.start()
    pb = ProgressBar(total=1)
    pb.start()

    # scan / stats loop {{{3
    print >> sys.stderr, '\n'
    sigints = 0
    sigintt = time.time()
    while not at.finished and at.exception is None:
        time.sleep(1)
        stats = engine.stats()
        if not stats['records_parsed']:
            continue

        if args.progress:
            pb.update(stats['progress'])
            print >> sys.stderr, str(pb),

        if args.coverage:
            means = sorted([n/len(analyser[i])
                    for i, n in enumerate(stats['nseqbasehits'])])

            if means and means[len(means)/2] > args.coverage:
                print >> sys.stderr
                lo.info('aborting scanning: median of coverage %d > %d'%(
                        means[len(means)/2], args.coverage))
                engine.stop()
                break

        # <CTRL-C> : output additional information
        if stats['sigints'] > sigints:

            # 2nd time : cancel scanning
            if time.time() - sigintt < 2.:
                print >> sys.stderr, '\n\n*** caught multiple <CTRL-C> within 2s : abort scanning ***'
                engine.abort()
                at.join()
                break

            print
            print
            print TextHist(title='readlengths').draw(stats['readlengths'], indexed=True)

            print
            print 'mean coverages'
            print '--------------'
            means = sorted([n/len(analyser[i])
                    for i, n in enumerate(stats['nseqbasehits'])])
            print TextHist().draw(sorted(means), indexed=False)

            sigints = stats['sigints']
            sigintt = time.time()

    at.join()
    if at.exception:
        lo.error('could not scan %s : %s [%s]'%(args.fastq, str(at.exception), at.traceback))
        sys.exit(ERROR_FASTQ_FORMAT_ERROR)

    print >> sys.stderr
    lo.info('performed scanning of %.2f%% in %.3f seconds'% (
            1e2*stats['progress'], time.time()-t0))

    # save to file {{{2
    analyser.update_testsuites()

    data = analyser.encode(hits=args.hits)
    j = codecs.open(args.json, 'w', 'utf-8')
    json_dump(data, j)


# show {{{1

def show(args):

    fastq = Fastq(args.file)

    if args.quality:
        Amin = fastq.Q2A(args.quality)
        n = args.number
        points = args.points
        lo.info('determining readlengths with quality>=%d of %s '
                'by reading %d records at %d points'%(
                args.quality, args.file, n, points))
        rls = fastq.lengths(Amin, n=n, points=points)

        hist = TextHist()
        print hist.draw(sorted(rls))

    if args.info:
        print 'dQ=' + str(fastq.dQ)
        print 'variants=' + str(fastq.variants)
        print 'readlength=' + str(fastq.readlength)
        print 'records_approx=' + str(fastq.records_approx)


# update {{{1

def update(args, testsuites):

    if args.fastq:
        lo.warning('re-reading of hits not currently implemented')

    f = open(args.json)
    data = json.load(f)

    analyser = analyse.Analyser()
    analyser.decode(testsuites, data)
    analyser.update_testsuites()

    # save results back to .json
    data = analyser.encode(hits = analyser.hits is not None)
    j = codecs.open(args.json, 'w', 'utf-8')
    json.dump(data, j, indent=2)


# illustrate {{{1

def illustrate(args, testsuites):

    analyser = analyse.Analyser()
    lo.info('loading json-file args.file')
    analyser.decode(testsuites, json.load(file(args.file)))
    lo.info('updating testsuites')
    analyser.update_testsuites()

    if args.readlengths:
        rls = analyser.stats['readlengths']

        hist = TextHist()
        print hist.draw(rls, indexed=True)

    if args.coverage:
        for name, testsuite in analyser.testsuites.items():
            print name + ':'
            for test in testsuite.tests:
                print '  - %s : %s' % (test, analyser[test])
            print

    if args.results:
        for testsuite, results in analyser.results.items():
            print '\n'+testsuite
            print '-'*len(testsuite)
            pprint(results)


# version {{{1

def version(args, testsuites):
    print VERSION


# gui {{{1

def gui(args, testsuites):
    # only import Tkinter etc now
    import Tkinter
    from kvarq.gui.main import MainGUI
    MainGUI(testsuites=testsuites)
    Tkinter.mainloop()



# info {{{1

def info(args, testsuites):
    print 'version=' + VERSION
    print 'testsuites=' + ','.join([
            '%s-%s' % (name, testsuite.version)
            for name, testsuite in testsuites.items()])
    print 'sys.prefix=' + sys.prefix


# explorer {{{1

def explorer(args, testsuites):
    import Tkinter as tk
    from kvarq.gui.explorer import DirectoryExplorer, JsonExplorer
    if os.path.isdir(args.explorable):
        DirectoryExplorer(args.explorable, testsuites=testsuites)
    else:
        JsonExplorer(args.explorable, testsuites=testsuites)
    tk.mainloop()


# parser {{{1

parser = argparse.ArgumentParser(description='''

        analyse .fastq file and report specific mutations in a .json file;
        additional output is displayed on stdout and log information is printed
        on stderr

    ''')

subparsers = parser.add_subparsers(help='main command to execute')

parser.add_argument('-d', '--debug', action='store_true',
        help='output log information at a debug level')
parser.add_argument('-q', '--quiet', action='store_true',
        help='only output warnings/errors to stderr/log')
parser.add_argument('-x', '--excepthook', action='store_true',
        help='catch exception and launch debugger')
parser.add_argument('-l', '--log',
        help='append log to specified file (similar to redirecting stderr, but without progress bar)')
parser.add_argument('-t', '--testsuites', action='append',
        help='load additional testsuites from specified files/directory (can be specified several times). all python files in specified directory not beginning with `_\' are loaded. filenames can contain version information after a `-\'. testsuites replace testsuites with same name specified in preceeding command line arguments. COMPULSORY switch for some commands (scan, update, illustrate, explorer)')

# version {{{2
parser_version = subparsers.add_parser('version',
        help='show version info')
parser_version.set_defaults(func=version)

# scan {{{2
parser_scan = subparsers.add_parser('scan')
parser_scan.set_defaults(func=scan)

# output console
parser_scan.add_argument('-p', '--progress', action='store_true',
        help='shows progress bar on stdout while scanning')

# general
parser_scan.add_argument('-S', '--no-scan', action='store_true',
        help='instead of scanning the original file, the provided .json file from a previous scan result is used and the .json structures are re-calculated from the hit-list (see -h); can only be done if version matches and same test suit is used (see -g)')

# scan config
parser_scan.add_argument('-t', '--threads', action='store', type=int,
        default=default_config['threads'],
        help='number of threads for concurrent scanning (default: %d)' % default_config['threads'])
parser_scan.add_argument('-Q', '--quality', action='store', type=int,
        default=default_config['quality'],
        help='discard nucleotides with Q score inferior to this value (default=%d; i.e. p=0.05)' % default_config['quality'])
parser_scan.add_argument('-e', '--errors', action='store', type=int,
        default=default_config['errors'],
        help='maximal number of consecutive errors allowed when comparing base sequences (default=%d)' % default_config['errors'])
parser_scan.add_argument('-r', '--readlength', action='store', type=int,
        default=default_config['minimum readlength'],
        help='minimum read length (default=%d)' % default_config['minimum readlength'])
parser_scan.add_argument('-o', '--overlap', action='store', type=int,
        default=default_config['minimum overlap'],
        help='minimum read overlap (default=%d)' % default_config['minimum overlap'])
parser_scan.add_argument('-c', '--coverage', type=int,
        default=default_config['stop median coverage'],
        help='stop scanning when median coverage (including margins) is above specified value (default=%d) -- specify 0 to force scanning of entire file' % default_config['stop median coverage'])
parser_scan.add_argument('-1', '--no-reverse', action='store_true',
        help='do not scan for hits in reverse strand')

# output .json
parser_scan.add_argument('-f', '--force', action='store_true',
        help='overwrite any existing .json file')
parser_scan.add_argument('-H', '--hits', action='store_true',
        help='saves all hits in .json file; this way scan result can be re-used without (see -n)')

# main arguments
parser_scan.add_argument('fastq',
        help='name of .fastq file to scan')
parser_scan.add_argument('json',
        help='name of .json file to where results are stored (or loaded, see -S)')

# update {{{2
parser_update = subparsers.add_parser('update',
        help='update (re-calculate) testsuites based on coverages saved in .json file; result is stored in same file')
parser_update.set_defaults(func=update)

parser_update.add_argument('json',
        help='name of .json file to update')
parser_update.add_argument('fastq', nargs='?',
        help='also re-calculate coverages with .fastq file specified (when .fastq file is not specified, coverages are taken from .json)')

# show {{{2
parser_show = subparsers.add_parser('show',
        help='show some information about a .fastq file')
parser_show.set_defaults(func=show)

parser_show.add_argument('-n', '--number', action='store', default=10000, type=int,
        help='number of records to read (applies to -Q)')
parser_show.add_argument('-p', '--points', action='store', default=10, type=int,
        help='number of points in file where to sample (spaced evenly; applies to -Q)')

parser_show.add_argument('-Q', '--quality', action='store', default=0, type=int,
        help='show histogram of readlengths with given quality cutoff (see also -n, -o)')
parser_show.add_argument('-i', '--info', action='store_true',
        help='output some information about FastQ file')

parser_show.add_argument('file',
        help='name of .fastq file to analyze')


# illustrate {{{2
parser_illustrate = subparsers.add_parser('illustrate',
        help='illustrate some information contained in a .json file (previously generated using the "scan" command)')
parser_illustrate.set_defaults(func=illustrate)

parser_illustrate.add_argument('-l', '--readlengths', action='store_true',
        help='show a histogram of readlengths')
parser_illustrate.add_argument('-c', '--coverage', action='store_true',
        help='show tests/coverages sorted by testsuite')

parser_illustrate.add_argument('-r', '--results', action='count',
        help='shows results of analyses')

parser_illustrate.add_argument('file',
        help='name of .json file to illustrate')


# gui {{{2
parser_gui = subparsers.add_parser('gui',
        help='start GUI')
parser_gui.set_defaults(func=gui)

# info {{{2
parser_info = subparsers.add_parser('info',
        help='show infos about kvarq')
parser_info.set_defaults(func=info)

# explorer {{{2
parser_explorer = subparsers.add_parser('explorer',
        help='launches the directory/json explorer')
parser_explorer.add_argument('explorable',
        help='directory/.json file to explore')
parser_explorer.set_defaults(func=explorer)


# __main__ {{{1

def main():
    args = parser.parse_args(sys.argv[1:])

    assert not (args.debug and args.quiet), \
            'make up your mind: debug OR normal OR quiet'

    if args.debug:
        set_debug()
    if args.quiet:
        set_warning()
    if args.log:
        appendlog(args.log)
    if args.excepthook:
        sys.excepthook = traceit

    testsuites = dict()

    if args.testsuites:
        for path in args.testsuites:
            if os.path.isdir(path):
                fnames = [fname for fname in glob.glob(os.path.join(path, '*.py'))
                        if not os.path.basename(fname)[0] == '_']
            else:
                fnames = [path]

            for fname in fnames:
                try:
                    name, testsuite = genes.load_testsuite(fname)
                    if name in testsuites:
                        lo.info('replaced testsuite "%s" v%s -> v%s from "%s"' %
                                (name, testsuites[name].version, testsuite.version, fname))
                    else:
                        lo.info('loaded testsuite "%s" from "%s"' % (name, fname))

                    testsuites[name] = testsuite

                except genes.TestsuiteLoadingException, e:
                    lo.error('could not load testsuite from "%s" : %s' % (fname, e))

    if args.func in [scan, update, illustrate, explorer]:
        if not testsuites:
            invalid_args('must specify at least one testsuite')
        args.func(args, testsuites=testsuites)
    elif args.func in [info, gui]:
        args.func(args, testsuites=testsuites)
    else:
        args.func(args)

if __name__ == "__main__":
    main()

