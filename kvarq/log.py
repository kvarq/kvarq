
from kvarq.util import is_exe_gui, is_app

import logging
import sys
import traceback
import os.path
import time
import functools
import re

class ColoredFormatter(logging.Formatter):

    def format(self, record):
        ret = super(ColoredFormatter, self).format(record)

        reset = '\033[m'
        bold = '\033[1m'
        ok = '\033[92m' # fg=green
        error = '\033[97;101m' # fg=white bg=red

        m = re.match('(\\[.*?\\] )(-INFO-)( .*)', ret)
        if m:
            return m.group(1) + ok + m.group(2) + reset + m.group(3)
        m = re.match('(\\[.*?\\] )(-WARNING-|-ERROR-)( .*)', ret)
        if m:
            return m.group(1) + error + m.group(2) + reset + bold + m.group(3) + reset

        return ret

lo = logging.getLogger('kvarq')
ft = logging.Formatter('[%(asctime)s] -%(levelname)s- %(filename)s:%(lineno)d(%(funcName)s) :: %(message)s')
cft = ColoredFormatter('[%(asctime)s] -%(levelname)s- %(filename)s:%(lineno)d(%(funcName)s) :: %(message)s')

def color_wrap(handler, emit_bak):
    def emit(record):
        import pprint
        pprint.pprint(record)
    return emit

logfn = None # use this to check whether log goes to file
if is_exe_gui() or is_app():
    # log to file in home directory if GUI is started from .exe / .app
    logfn = os.path.join(os.path.expanduser('~'), 'kvarq.log')
    ch = logging.FileHandler(logfn)
    ch.setFormatter(ft)
else:
    # log to stderr if started from command line
    ch = logging.StreamHandler(sys.stderr)
    if sys.platform != 'win32' and hasattr(sys.stderr, 'isatty') and sys.stderr.isatty():
        ch.setFormatter(cft)
    else:
        ch.setFormatter(ft)

lo.setLevel(logging.INFO)
ch.setLevel(logging.DEBUG)
lo.addHandler(ch)

def set_debug():
    lo.setLevel(logging.DEBUG)
def set_info():
    lo.setLevel(logging.INFO)
def set_warning():
    lo.setLevel(logging.WARNING)

def appendlog(fname):
    fh = logging.FileHandler(fname)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(ft)
    lo.addHandler(fh)

tictocs = {}
def tic(name):
    tictocs.setdefault(name, []).append([time.time()])
def toc(name):
    l = tictocs[name][-1]
    l.append(time.time())
    lo.debug('toc-tic %s : %.2f ms'%(name, 1e3*(l[1]-l[0])))

def tictoc(name):
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            tic(name)
            ret = f(*args, **kwargs)
            toc(name)
            return ret
        return wrapper
    return decorator

def format_traceback(exc_info):
    return ' -> '.join([
        '%s:%d' % (frame[0], frame[1])
        for frame in traceback.extract_tb(exc_info[2])])

