
from kvarq.log import lo, tic, toc
from kvarq import VERSION, genes, engine, analyse
from kvarq.fastq import Fastq, FastqFileFormatException
from kvarq.gui.explorer import JsonExplorer
from kvarq.engine import Hit
from kvarq.util import ProgressBar
from kvarq.gui.util import open_help, ThemedTk, askopenfilename

import Tkinter as tk
import tkFont
import tkFileDialog
import tkMessageBox

import os.path
import glob
import sys
import logging
import threading
import time
import json


class AnalyseThread(threading.Thread):

    def __init__(self, analyser, fastq, testsuites):
        super(AnalyseThread, self).__init__(name='analyse-thread')
        self.analyser = analyser
        self.fastq = fastq
        self.testsuites = testsuites
        self.finished = False
        self.exception = None
        self.stopped = False

    def run(self):
        try:
            self.analyser.scan(self.fastq, self.testsuites)
            self.finished = True
        except Exception, e:
            self.exception = e

    def stop(self):
        engine.stop()
        self.stopped = True


class SimpleGUI(ThemedTk):

    def __init__(self, settings):
        self.settings = settings

        self.fastqi = -1
        self.analysers = {}
        self.fastqs = self.askfastqs()
        if self.fastqs is not None:
            self.init_gui()
            self.next_fastq()


    def init_gui(self):
        ThemedTk.__init__(self, title='scan .fastq files')

        self.bind('<Destroy>', self.destroy_cb)

        self.analyser = None
        self.running = False
        self.save_hits = False

#        self.menu = tk.Menu(self)
#        filemenu = tk.Menu(self.menu)
#        self.menu.add_cascade(label='KvarQ', menu=filemenu)
#        filemenu.add_command(label='Help', command=open_help)
#        self.config(menu=self.menu)

        frame = tk.Frame(self)
        self.fname = tk.Label(frame)
        self.fname.pack(side=tk.LEFT)
        frame.pack(side=tk.TOP, expand=False, fill=tk.X)

        frame = tk.Frame(self)
        self.running = False
        self.start = tk.Button(frame, text='start', command=self.startstop, state=tk.DISABLED)
        self.start.pack(side=tk.LEFT)
        self.pb = ProgressBar(total=1., r='')
        self.pb_longest = 0
        self.pblabel = tk.Label(frame, text=' ' * 70, font=self.monospace)
        self.pblabel.pack(side=tk.LEFT)
        frame.pack(side=tk.TOP, expand=False, fill=tk.X)

        frame = tk.Frame(self)
        self.show = tk.Button(frame, text='show', command=self.show_cb, state=tk.DISABLED)
        self.show.pack(side=tk.LEFT)
        self.save = tk.Button(frame, text='save', command=self.save_cb, state=tk.DISABLED)
        self.save.pack(side=tk.LEFT)
        frame.pack(side=tk.TOP, expand=False, fill=tk.X)

        if len(self.fastqs) > 1:
            self.show.config(text='show last')
            self.save.config(text='save all')

        self.resizable(0, 0)


    def askfastqs(self):

        ret = None
        while ret is None:

            fastqs = askopenfilename(
                    multiple=True,
                    filetypes=[
                        ('fastq files', '*.fastq'),
                        ('compressed fastq files', '*.fastq.gz')],
                    title='select .fastq files to analyze')

            if not fastqs:
                return None

            ret = fastqs
            for fastq in fastqs:
                if not os.path.isfile(fastq) and fastq.lower().endswith('.fastq'):
                    tkMessageBox.showerror(
                            'invalid file selected',
                            '"%s" does not seem to be a .fastq file' % fastq)
                    ret = None
                    break

        return ret


    def has_more_fastq(self):
        return self.fastqi + 1 < len(self.fastqs)

    def next_fastq(self):
        n = len(self.fastqs)

        while True:
            self.fastqi += 1
            if self.fastqi >= n:
                return False

            try:
                self.fastq = Fastq(self.fastqs[self.fastqi])
            except FastqFileFormatException, e:
                lo.error('cannot load file %s : %s'%(self.fastqs[self.fastqi], str(e)))
                if n == 1:
                    tkMessageBox.showerror(
                            'invalid .fastq file',
                            'the selected file cannot be parsed : ' + str(e))
                continue

            if n == 1:
                self.fname.config(text=self.fastq.fname)
            else:
                self.fname.config(text=self.fastq.fname + ' (file %d/%d)'%(self.fastqi+1, n))

            self.start.config(state=tk.NORMAL)

            return True


    def startstop(self):
        if not self.running:
            # "start"
            self.analyser = analyse.Analyser()

            engine.config(
                    nthreads=self.settings.config['threads'],
                    maxerrors=self.settings.config['errors'],
                    minreadlength=self.settings.config['minimum readlength'],
                    minoverlap=self.settings.config['minimum overlap'],
                    Amin=self.fastq.Q2A(self.settings.config['quality']),
                    Azero=self.fastq.Azero
                )

            self.at = AnalyseThread(self.analyser, self.fastq, self.settings.get_testsuites())
            self.t0 = time.time()
            self.at.start()
            self.pb.start()
            self.after(100, self.update)
            lo.info('start scanning %s (%d MB)'%(
                    self.fastq.fname,
                    os.path.getsize(self.fastq.fname)/1024**2))
            self.running = True

            self.start.config(text='stop')

        else:
            # "stop"
            if self.fastqi+1 < len(self.fastqs):
                if self.next_fastq():
                    self.start.config(text='start next')
                else:
                    self.start.config(state=tk.DISABLED)
            else:
                self.start.config(state=tk.DISABLED)
            self.running = False

    def finish_scanning(self):
        lo.info('analyzing data...')
        self.analyser.update_testsuites()
        lo.info('done analyzing data')
        self.analysers[self.fastqs[self.fastqi]] = self.analyser
        self.show.config(state=tk.NORMAL)
        self.save.config(state=tk.NORMAL)

    def update(self):

        if not self.running:
            self.at.stop()
            self.at.join()
            lo.info('STOPPED scanning via GUI after %.3f seconds'% (time.time()-self.t0))
            self.finish_scanning()
            self.running = False
            self.at = None
            return

        stats = engine.stats()
        self.pb.update(stats['progress'])
        pb_str = str(self.pb)
        self.pb_longest = max(self.pb_longest, len(pb_str)) # prevent too much window resizing
        self.pblabel.config(text=('{:<%d}' % self.pb_longest).format(pb_str))

#        if self.settings.config['stop median coverage'] and not self.at.stopped:
#            means = sorted([n/len(self.analyser[i])
#                    for i, n in enumerate(stats['nseqbasehits'])])
#
#            if means and means[len(means)/2] > self.settings.config['stop median coverage']:
#                lo.info('aborting scanning: median of coverage %d > %d'%(
#                        means[len(means)/2], self.settings.config['stop median coverage']))
#                self.at.stop()

        if self.at.finished or self.at.exception:
            self.at.join()
            self.start.config(state=tk.DISABLED)
            if self.at.finished:
                lo.info('finished scanning after %.3f seconds'% (time.time()-self.t0))
                self.pblabel.config(text=str(self.pb)[:str(self.pb).index(']')+1]+' -- done')
                self.finish_scanning()
            if self.at.exception:
                lo.error('could not scan %s : %s'%(self.fastq.fname, str(self.at.exception)))
            self.running = False
            self.at = None
            if self.next_fastq():
                self.startstop()
            return

        self.after(100, self.update)

    def show_cb(self):
        if self.analyser.results is None:
            tkMessageBox.showinfo('no results yet', 'please stop/finish the scanning first')
            return
        explorer = JsonExplorer(self.analyser, testsuites=self.settings.get_testsuites())

    def save_cb(self):

        if len(self.analysers) == 1:
            jf = tkFileDialog.asksaveasfile(
                    parent=self,
                    initialfile=os.path.splitext(os.path.basename(
                            self.fastq.fname))[0] + '.json',
                    initialdir=os.path.dirname(self.fastq.fname),
                    defaultextension='.json', 
                    title='select .json to store results of scan')
            if not jf:
                return
            tic('dumping json')
            data = self.analyser.encode(hits=self.save_hits)
            json.dump(data, jf, indent=2)
            toc('dumping json')

        else:
            jd = tkFileDialog.askdirectory(
                    parent=self,
                    title='select directory to store .json files')
            if not jd:
                return
            for fastq, analyser in self.analysers.items():
                base = os.path.splitext(os.path.basename(fastq))[0]
                while True:
                    jsonfn = os.path.join(jd, base + '.json')
                    if not os.path.exists(jsonfn):
                        break
                    base += '_'
                lo.info('saving to ' + jsonfn)
                tic('dumping json')
                data = analyser.encode(hits=self.save_hits)
                json.dump(data, file(jsonfn, 'w'), indent=2)
                toc('dumping json')

    def destroy_cb(self, x=None):
        if self.running:
            self.at.stop()
            self.at.join()
            lo.removeHandler(self.log_handler)


if __name__=='__main__':
            
    from kvarq.gui.settings import Settings
    win = SimpleGUI(Settings())
    tk.mainloop()

