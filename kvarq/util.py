
import time, sys
from math import log
import threading
import os.path
import json
import re


class ProgressBar(object):

    ''' [==========>   32%               ] 1m 30s / 4m 14s
        [============= 65% =====>        ] 3m 02s / 4m 14s '''

    def __init__(self, total, width=50, ETA=True, r=None, ETAbuflen=100):
        self.total = total
        self.started = None
        self.width = width
        self.ETA = ETA
        self.ETAbuflen = ETAbuflen
        self.ETAbuf = []

        if r is None:
            if sys.platform.lower().startswith('win'):
                self.r = '\r'
            else:
                self.r = '\n' + chr(27) + '[A'
        else:
            self.r = r

    def start(self):
        self.started = time.time()

    def update(self, done):
        self.ETAbuf.append((done, time.time()))
        while len(self.ETAbuf)>self.ETAbuflen:
            del self.ETAbuf[0]

    def fmt_secs(self, secs):
        if secs > 3600:
            return '%dh %dm %ds'%(
                    int(secs/3600),
                    int((secs%3600)/60),
                    int(secs%60))
        elif secs>60:
            return '%dm %ds'%(
                    int(secs/60),
                    int(secs%60))
        else:
            return '%ds'%int(secs)

    def get_ETA(self):
        if self.ETA:
            secs_done = self.when - self.started
            if self.ETAbuf and secs_done>10:
                secs_total = sum([self.total/done*(when-self.started)
                        for done, when in self.ETAbuf if done]) / len(self.ETAbuf)
                if secs_total>120:
                    secs_total = int(secs_total/10) * 10
                return ' ' + self.fmt_secs(secs_done) + ' / ' + \
                        self.fmt_secs(secs_total) + '   '
            else:
                return ' ' + self.fmt_secs(secs_done) + ' / ???'
        else:
            return ''

    def barize(self, p, maxbars):
        return '='*int(p*maxbars) + '>' + ' '*int((1-p)*maxbars)

    def __str__(self):
        if not self.started:
            self.start()
        secs_done = self.started - self.when
        p = min(0.9999, float(self.done) / self.total)
        percents = '%3d%% ' % int(100*p)
        maxbars = self.width - len(percents) - len('[]')

        if p<0.5:
            return self.r + '[' + self.barize(2*p, maxbars/2) + percents + \
                    ' '*int(maxbars/2) + ']' + self.get_ETA()
        else:
            return self.r + '[' + '='*int(maxbars/2) + percents + \
                     self.barize(2*(p-0.5), maxbars/2) + ']' + self.get_ETA()

    @property
    def done(self):
        if self.ETAbuf: return self.ETAbuf[-1][0]
        return None
    @property
    def when(self):
        if self.ETAbuf: return self.ETAbuf[-1][1]
        return None

    @classmethod
    def run_watched(cls, f_cb, progress_cb, dt=1.):
        ''' runs f_cb and shows progress bar based on progress_cb()
            in another thread '''

        class ProgressThread(threading.Thread):

            def __init__(self, f_cb):
                super(ProgressThread, self).__init__(name='progressbar-thread')
                self.pb = ProgressBar(total=1.)
                self.ret = None
                self.f_cb = f_cb
                self.done = False

            def run(self):
                self.ret = self.f_cb()
                self.done = True

        pt = ProgressThread(f_cb)
        pb = cls(total=1.)
        pt.start()
        pb.start()
        print >> sys.stderr

        while not pt.done:
            pb.update(progress_cb())
            print >> sys.stderr, str(pb)
            time.sleep(dt)

        return pt.ret




class TextHist:

    ''' outputs a text histogram in the form of::

            [  0-100]  33(10%)********
            [100-200] 185(50%)****************************************
            [200-300]  66(20%)****************
            [300-400]  40(15%)************
            [400-500]  40(15%)************
            [500-600]   0( 0%)

            totaling 330, average 172.35 '''

    def __init__(self, bins=15, width=65, title=None):
        self.bins = bins
        self.width = width
        self.title = title


    def draw(self, data, indexed=False):
        ''' data is already sorted

              - indexed=``False`` : data is array containing actual values
              - indexed=``True`` : each entry in data represents number of
                occurences of its index '''

        if not data:
            return 'no data --> CANNOT GENERATE HISTOGRAM'
        if indexed:
            bw = len(data)/float(self.bins)
            N = int(log(len(data))/log(10)) +1
        else:
            if data[-1]==0:
                return '*** cannot draw histogram because all data zero ***'
            bw = (data[-1]-data[0])/float(self.bins)
            N = int(log(data[-1])/log(10)) +1
        if not bw:
            return 'bw=0 --> CANNOT GENERATE HISTOGRAM'
        n = int(log(bw)/log(10)) -1
        bw = int(bw/10**n) * 10**n
        bw = max(bw, 1.)

        xs = []
        i = bi = x = mx = sx = s = 0
        while i<len(data):
            if not indexed and (data[i]>(bi+1)*bw) or \
                    indexed and (i>(bi+1)*bw):
                xs.append(x)
                sx += x
                if x > mx: mx = x
                x = 0
                bi += 1
            else:
                if indexed:
                    x += data[i]
                    s += data[i]*i
                else:
                    x += 1
                    s += data[i]
                i += 1
        if x:
            xs.append(x)

        fmt = '[%%%dd-%%%dd] %%%dd (%%2d%%%%)'%(
                max(N,4), max(N,4), int(log(mx)/log(10)) +1)
        ret = ''
        if self.title:
            ret += self.title + '\n' + '-'*(len(self.title)) + '\n'
        for bi, x in enumerate(xs):
            ret += fmt%(bi*bw, (bi+1)*bw, x, int(100*x/sx))
            ret += '*'*int(self.width*x/mx) + '\n'

        if indexed:
            avg = float(s)/sum(data)
        else:
            avg = float(s)/len(data)
        ret += 'totaling %d, average %.2f'%(mx, avg)

        return ret

def get_root_path(*parts):
    if is_exe():
        return os.path.join(sys.prefix, *parts)
    elif is_app():
        return os.path.join(sys.prefix, os.path.pardir, os.path.pardir,
                os.path.pardir, *parts)
    else:
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
        return os.path.join(root, *parts)

def is_exe_console():
    return hasattr(sys, 'frozen') and sys.frozen=='console_exe' # pylint: disable=E1101

def is_exe_gui():
    return hasattr(sys, 'frozen') and sys.frozen=='windows_exe' # pylint: disable=E1101

def is_exe():
    return is_exe_gui() or is_exe_console()

def is_app():
    return hasattr(sys, 'frozen') and (
            sys.frozen=='macosx_app') # pylint: disable=E1101


def json_dump(data, fd, indent=2, max_indent_level=2):
    ''' 
    :param data: python object to dump
    :param fd: open file for writing data to
    :param indent: number of spaces per level of indentation
    :param max_indent_level: number of levels to indent

    writes data "nicely formatted" to specified file
    '''

    ii = indent * max_indent_level
    re1 = re.compile('^([\\[{,]? ?)\n {%d,}' % (ii + 1), re.MULTILINE)
    re2 = re.compile('^\n {%d}$' % ii, re.MULTILINE)
    spacer = ''

    for chunk in json.JSONEncoder(indent=2).iterencode(data):
        if re2.match(chunk):
            spacer = chunk
        else:
            if spacer and chunk not in ['}', ']', ',']:
                fd.write(spacer)
            spacer = ''
            fd.write(re1.sub('\\1', chunk))


class csv_xls_writer:

    @classmethod
    def add_extension(cls, fname):
        try:
            import xlwt
            return os.path.splitext(fname)[0] + '.xls'
        except ImportError:
            return os.path.splitext(fname)[0] + '.csv'

    def __init__(self, fname, autoflush=True, sheet_name='exported data'):

        ''' create a new ``csv_xls_writer`` instance

            if you specify a file name that ends with ``.xls``, it will be
            renamed into ``.csv`` if the module ``xlwt`` is not installed. '''

        self.fname = fname

        if fname.endswith('.csv'):
            import csv
            self.csv = csv.writer(file(self.fname, 'w'))
            self.xls = None

        elif fname.endswith('.xls'):
            try:
                import xlwt
            except ImportError:
                import csv
                self.fname = self.fname[:-4] + '.csv'
                self.csv = csv.writer(file(self.fname, 'w'))
                self.xls = None
                return
            self.easyxf = xlwt.easyxf
            self.wb = xlwt.Workbook()
            self.ws = self.wb.add_sheet(sheet_name)
            self.row = 0
            self.autoflush = autoflush
            self.csv = None

        else:
            raise IOError('can only export data to .csv or .xls')

    def writerow(self, row, colors=None):
        if self.csv:
            self.csv.writerow(row)
        else:
            for col, value in enumerate(row):
                # excel colors : http://dmcritchie.mvps.org/excel/colors.htm
                if colors and col in colors:
                    st = self.easyxf('pattern: pattern solid;')
                    st.pattern.pattern_fore_colour = colors[col]
                    self.ws.write(self.row, col, value, st)
                else:
                    self.ws.write(self.row, col, value)
            self.row += 1
            if self.autoflush:
                self.flush()

    def flush(self):
        if self.csv: return
        self.wb.save(self.fname)

