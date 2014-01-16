
from collections import Counter
from math import log
import Tkinter as tk
import tkFont
import tkMessageBox

from kvarq.gui.util import ThemedTk


class IndexableDisplay(object):

    def __init__(self, parent, data):
        self.frame = tk.Frame(parent)
        self.grid = self.frame.grid
        self.pack = self.frame.pack

        self.canvas = tk.Canvas(self.frame)
        self.canvas.pack(fill='both', expand=1)
        self.canvas.bind('<Configure>', self.update)
        self.canvas.bind('<Motion>', self.motion)

        self.label = tk.Label(self.frame)
        self.label.pack()

        self.margin = [20, 20, 20, 20]
        self.config = {
                'ticks' : {
                    'length' : 8,
                    'width' : 2,
                    'font' : tkFont.Font(self.canvas, size=10),
                },
            }

        self.data = data

    def is_inside(self, x, y):
        return (x >= self.margin[3] and
                x < self.canvas.winfo_width() - self.margin[1] and
                y >= self.margin[0] and
                y < self.canvas.winfo_height() - self.margin[2])

    def get_index(self, x):
        return len(self.data) * (x - self.margin[3]) / self.width()

    def width(self):
        return self.canvas.winfo_width() - self.margin[1] - self.margin[3]
    def height(self):
        return self.canvas.winfo_height() - self.margin[0] - self.margin[2]

    def create_rectangle(self, x1, y1, w, h, *a, **b):
        self.canvas.create_rectangle(self.margin[3] + x1,
                self.canvas.winfo_height() - self.margin[2] - y1,
                self.margin[3] + x1 + w,
                self.canvas.winfo_height() - self.margin[2] - y1 - h,
                *a, **b)

    def create_line(self, x1, y1, x2, y2, *a, **b):
        self.canvas.create_line(self.margin[3] + x1,
                self.canvas.winfo_height() - self.margin[2] - y1,
                self.margin[3] + x2,
                self.canvas.winfo_height() - self.margin[2] - y2,
                *a, **b)

    def xticks(self, ticks, labels=None):
        if not labels:
            labels = [str(label) for label in ticks]

        w = self.width()
        ch = self.canvas.winfo_height()
        for i, tick in enumerate(ticks):
            x = int(w * (tick + 0.5) / len(self.data))

            self.create_line(x, 0, x, self.config['ticks']['length'],
                    fill='black', width=self.config['ticks']['width'])
            self.canvas.create_text(self.margin[3] + x, ch - self.margin[2] / 2,
                    text=labels[i])

    def transform_ylabel(self, number):
        if number > 1e9:
            return '%.2fG' % (number / 1e9)
        if number > 1e6:
            return '%.2fM' % (number / 1e6)
        if number > 1e3:
            return '%.2fk' % (number / 1e3)
        return str(number)

    def yticks(self, ticks, labels=None):
        if not labels:
            labels = [self.transform_ylabel(label) for label in ticks]

        h = self.height()
        ch = self.canvas.winfo_height()
        for i, tick in enumerate(ticks):
            y = h * tick / max(self.data)

            self.create_line(0, y, self.config['ticks']['length'], y,
                    fill='black', width=self.config['ticks']['width'])
            self.canvas.create_text(self.margin[3] / 2, ch - self.margin[2] - y,
                    text=labels[i])

    def update(self):
        pass

    def motion(self, e=None):
        pass


class CoverageDisplay(IndexableDisplay):

    def __init__(self, parent, pos0, coverage):
        super(CoverageDisplay, self).__init__(parent, coverage.coverage)

        self.colors = {
                'coverage' : 'gray',
                'bases' : {
                    'A' : 'blue',
                    'C' : 'cyan',
                    'G' : 'green',
                    'T' : 'red',
                },
            }

        self.coverage = coverage
        self.pos0 = pos0

    def motion(self, e):

        if not self.is_inside(e.x, e.y):
            self.label.config(text='(mouse outside graph)')
            return

        i = self.get_index(e.x)

        if self.pos0:
            text = 'pos=' + str(self.pos0 + i - self.coverage.start)
        else:
            text = 'pos=?'
        text += ' coverage=' + str(self.coverage.coverage[i])
        mutations = self.coverage.mutations.get(i)
        if mutations:
            text += ' mutations=' + ', '.join(['%dx %s' % (v, k)
                    for k, v in Counter(mutations).items()])

        self.label.config(text=text)

    def update(self, e=None):

        self.canvas.delete('all')
        w, h = self.width(), self.height()

        x = self.data
        self.hitheight = min(3., float(h) / max(x))
        for i in range(len(x)):
            rx = w * i / len(x)
            rh = int(x[i] * self.hitheight)
            rw = w * (i + 1) / len(x) - w * i / len(x)
            self.create_rectangle(rx, 0, rw, rh,
                    fill=self.colors['coverage'], outline='')

            for j, b in enumerate(sorted(self.coverage.mutations.get(i, ''))):
                my = int(j * self.hitheight)
                mh = int((j + 1) * self.hitheight) - int(j * self.hitheight)
                self.create_rectangle(rx, my, rw, mh,
                        fill=self.colors['bases'][b], outline='')

        for i in (self.coverage.start, self.coverage.stop):
            lx = w * i / len(x)
            self.create_line(lx, 0, lx, h, fill='red', width=2)

        mean = self.coverage.mean(include_margins=False)
        std = self.coverage.std(include_margins=False)
        for j, dash in ((mean - std, (2, 5)), (mean, (5, 3)), (mean + std, (2, 5))):
            ly = int(j * self.hitheight)
            self.create_line(0, ly, w, ly, fill='black', dash=dash)

        self.create_rectangle(0, 0, w, h, fill='', outline='black')


class CoverageWindow:

    def __init__(self, test, coverage):

        if not coverage.coverage or not max(coverage.coverage):
            tkMessageBox.showinfo('cannot show coverage',
                    'this template has an empty coverage')
            return

        title = str(test)
        self.win = ThemedTk(title=title, esc_closes=True, geometry=(800, 600))

        pos0 = None
        if hasattr(test.template, 'start'):
            pos0 = test.template.start
        self.cd = CoverageDisplay(self.win, pos0, coverage)
        self.cd.pack(fill='both', expand=1)


class ReadlengthDisplay(IndexableDisplay):

    def __init__(self, parent, data):
        super(ReadlengthDisplay, self).__init__(parent, data)
        self.margin[3] = 60

    def motion(self, e):

        if not self.is_inside(e.x, e.y):
            self.label.config(text='(mouse outside graph)')
            return

        i = self.get_index(e.x)

        text = '%d reads with length=%d' % (self.data[i], i)

        self.label.config(text=text)

    def update(self, e=None):

        self.canvas.delete('all')
        w, h = self.width(), self.height()

        x = self.data
        first = last = greatest = None
        for i, rl in enumerate(x):
            if x[i] > 0 and first is None: first = i
            if x[i] > 0: last = i
            if greatest is None or x[i] > greatest: greatest = x[i]

            rx = w * i / len(x)
            rh = int(float(h) * x[i] / max(x))
            rw = w * (i + 1) / len(x) - w * i / len(x)
            self.create_rectangle(rx, 0, rw, rh,
                    fill='gray', outline='')

        #self.create_rectangle(0, 0, self.width(), self.height(), fill='', outline='black')
        self.create_line(0, 0, w, 0, fill='black')
        self.create_line(0, 0, 0, h, fill='black')

        self.draw_ticks(first, last, greatest)

    def draw_ticks(self, first, last, greatest):
        self.xticks([first, last, self.data.index(greatest)])
        self.yticks([greatest])


class ReadlengthWindow:

    def __init__(self, readlenghts):
        if readlenghts and max(readlenghts):
            title = 'Length of quality-cut reads'
            self.win = ThemedTk(title=title, esc_closes=True, geometry=(600, 400))
            self.cd = ReadlengthDisplay(self.win, readlenghts)
            self.cd.pack(fill='both', expand=1)
        else:
            tkMessageBox.showinfo('cannot show readlengths', 'no reads found')


class HitHistogramDisplay(IndexableDisplay):

    def __init__(self, parent, data, indexed=False, nbins=15):
        ''' - indexed=``False`` : data is array containing actual values
            - indexed=``True`` : each entry in data represents number of
              occurences of its index '''
        nbins += 1
        bins, binwidth = self.make_bins(data, indexed, nbins)

        super(HitHistogramDisplay, self).__init__(parent, bins)
        self.binwidth = binwidth
        self.smallest = min(data)
        self.largest = max(data)

    def make_bins(self, data, indexed, bins):
        
        if indexed:
            bw = len(data) / float(bins)
        else:
            data = sorted(data)
            bw = (data[-1] - data[0]) / float(bins)

        n = int(log(bw) / log(10)) - 1
        binwidth = int(bw / 10 ** n) * 10 ** n

        binwidth = max(1., binwidth)

        bins = []
        i = bi = x = mx = sx = s = 0
        while i<len(data):
            if not indexed and (data[i]>(bi+1)*binwidth) or \
                   indexed and (i>(bi+1)*binwidth):
                bins.append(x)
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
            bins.append(x)

        return bins, binwidth

    def motion(self, e):

        if not self.is_inside(e.x, e.y):
            self.label.config(text='(mouse outside graph)')
            return

        i = self.get_index(e.x)

        a, b = self.binwidth * i, self.binwidth * (i + 1) - 1
        if i == 0:
            a = self.smallest
        if i == len(self.data) - 1:
            b = self.largest

        text = '%d templates with %d' % (self.data[i], a)
        if a != b:
            text += '-' + str(b)
        text += ' hits'

        self.label.config(text=text)

    def update(self, e=None):

        self.canvas.delete('all')
        w, h = self.width(), self.height()

        x = self.data
        for i, rl in enumerate(x):
            rx = w * i / len(x)
            rh = int(float(h) * x[i] / max(x))
            rw = w * (i + 1) / len(x) - w * i / len(x)
            self.create_rectangle(rx, 0, rw, rh,
                    fill='gray', outline='')

        self.create_rectangle(0, 0, self.width(), self.height(), fill='', outline='black')


class HitHistogramWindow:

    def __init__(self, data, indexed=False):

        if data and max(data):
            title = 'Histogram of hits/template'
            self.win = ThemedTk(title=title, esc_closes=True, geometry=(600, 400))
            self.cd = HitHistogramDisplay(self.win, data, indexed)
            self.cd.pack(fill='both', expand=1)
        else:
            tkMessageBox.showinfo('cannot show hits/template', 'no hits to template found')


class MeanCoverageDisplay(HitHistogramDisplay):

    def motion(self, e):

        if not self.is_inside(e.x, e.y):
            self.label.config(text='(mouse outside graph)')
            return

        i = self.get_index(e.x)

        text = '%d templates with medium coverage %d-%d' % (
                self.data[i], self.binwidth * i, self.binwidth * (i + 1))

        self.label.config(text=text)


class MeanCoverageWindow:

    def __init__(self, data, indexed=False):
        if data and max(data):
            title = 'Mean coverage of templates'
            self.win = ThemedTk(title=title, esc_closes=True, geometry=(600, 400))
            self.cd = MeanCoverageDisplay(self.win, data, indexed)
            self.cd.pack(fill='both', expand=1)
        else:
            tkMessageBox.showinfo('cannot show mean coverage', 'no hits to template found')


class SpoligoDisplay(ReadlengthDisplay):

    def __init__(self, parent, data):
        super(SpoligoDisplay, self).__init__(parent, data)
        self.margin[1] = self.margin[3] = 0

    def motion(self, e):

        if not self.is_inside(e.x, e.y):
            self.label.config(text='(mouse outside graph)')
            return

        i = self.get_index(e.x)

        text = 'spoligo%d : %d hits' % (i, self.data[i])

        self.label.config(text=text)

    def draw_ticks(self, *a, **b):
        pass


class SpoligoWindow:

    def __init__(self, spoligos):
        if spoligos and max(spoligos):
            title = 'Hits/spoligo'
            self.win = ThemedTk(title=title, esc_closes=True, geometry=(600, 400))
            self.cd = SpoligoDisplay(self.win, spoligos)
            self.cd.pack(fill='both', expand=1)
        else:
            tkMessageBox.showinfo('cannot show hits/spoligo', 'no spoligos found')


if __name__ == '__main__':

    import json

    from kvarq.genes import ancestor
    from kvarq.analyse import Coverage

    spacing = 25
    start = 761082
    stop = 761162
    seq = ancestor.seq(start, stop, spacing, spacing)
    coverage = Coverage(seq)
    coverage.deserialize('20-21-22-24-24-25-29-29-30-30-31-32-33-33-33-33-33-36-38-40-41-41-41-41-42-41-44-43-43-44-44-44-45-43-43-44-44-46-47-44-44-43-44-45-44-43-43-43-41-40-40-35-34-33-34-33-34-34-33-32-32-30-29-30-28-28-26-25-25-25-22-21-22-23-22-22-23-23-21-22-23-22-21-21-21-20-19-20-18-18-18-19-18-19-18-18-18-18-17-17-16-14-14-13-13-13-9-9-9-9-9-9-9-8-8-8-7-7-7-7-5-5-5-5-5-5-4-4-4-4-4 37[A]-42[A]-53[G]-78[T]-109[A]')

    win = tk.Tk()
    cd = CoverageDisplay(win, start, coverage)
    cd.pack(fill='both', expand=1)
    win.wm_geometry('800x600')


    HitHistogramWindow([1,1,2,4,5,3,3,1], indexed=True)
    ReadlengthWindow([0,0,0,1,1,2,1,2,5,2,3,0,2,2,1,0,0])
    tk.mainloop()


