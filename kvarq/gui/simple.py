
from kvarq.log import lo, tic, toc
from kvarq import VERSION, genes, engine, analyse
from kvarq.fastq import Fastq, FastqFileFormatException
from kvarq.gui.explorer import JsonExplorer
from kvarq.engine import Hit
from kvarq.util import ProgressBar
from kvarq.gui.util import open_help, ThemedTk, askopenfilename
from kvarq.testsuites import load_testsuites
from kvarq.config import config_params

import Tkinter as tk
import tkFont
import tkFileDialog
import tkMessageBox

import os
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
#            lo.debug('AnalyseThread : start')
            self.analyser.scan(self.fastq, self.testsuites)
#            lo.debug('AnalyseThread : stop')
            self.finished = True
        except Exception, e:
#            lo.debug('AnalyseThread : exception')
            self.exception = e

    def stop(self):
        engine.stop()
        self.stopped = True


class TestsuiteSelector(ThemedTk):

    def __init__(self, testsuite_paths):
        ThemedTk.__init__(self)

        label = tk.Label(self, text='select testsuites for scan:')
        label.pack(anchor='w')

        self.values = {}
        self.buttons = []
        for name in sorted(testsuite_paths):
            self.values[name] = False

            def make_toggle(name):
                def toggle(e=None):
                    self.values[name] = not self.values[name]
                return toggle

            button = tk.Checkbutton(self, text=name, command=make_toggle(name))
            button.pack(anchor='w')
            self.buttons.append(button)

        self.disabled = False
        self.closed = False
        self.protocol('WM_DELETE_WINDOW', self.closing)

    def closing(self, e=None):
        self.closed = True
        self.destroy()

    def selection(self):
        return [
                name for name, value in self.values.items()
                if value
            ]

    def disable(self):
        if not self.closed:
            for button in self.buttons:
                button.config(state='disabled')
        self.disabled = True


class SimpleGUI(ThemedTk):

    def __init__(self, settings, testsuites, testsuite_paths):
        self.settings = settings
        self.testsuites = testsuites # all loaded testsuites
        self.testsuite_paths = testsuite_paths

        self.fastqi = -1
        self.analysers = {}
        self.fastqs, self.paireds = self.askfastqs()

        if self.fastqs is not None:
            self.init_gui()
            self.next_fastq()


    def init_gui(self):
        ThemedTk.__init__(self, title='scan .fastq files')

        self.bind('<Destroy>', self.destroy_cb)

        self.selector = TestsuiteSelector(self.testsuite_paths)
        self.selected_testsuites = {} # testsuites selected for scanning

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
        self.protocol('WM_DELETE_WINDOW', self.closing)


    def closing(self, e=None):
        if not self.selector.closed:
            self.selector.destroy()
        self.destroy()


    def askfastqs(self):
        ''' asks user to select '*.fastq' and '*.fastq.gz' files and returns
            a list of fastq files and a list of booleans indicating whether
            the corresponding file should be scanned in 'paired' mode 
            (or ``None, None`` if canceled)'''

        fastqs = askopenfilename(
                initialdir=os.getcwd(),
                multiple=True,
                filetypes=[
                    ('fastq files', '*.fastq'),
                    ('compressed fastq files', '*.fastq.gz')],
                title='select .fastq files to analyze')

        if not fastqs:
            return None, None

        fastqs = sorted(fastqs)
        paireds = []
        idx = 0
        while idx < len(fastqs) - 1:
            base1 = fastqs[idx]  [:fastqs[idx]  .rindex('.fastq')]
            base2 = fastqs[idx+1][:fastqs[idx+1].rindex('.fastq')]
            if len(base1) > 2 and len(base2) > 2 and \
                    base1[-2:] == '_1' and base2[-2:] == '_2' and \
                    base1[:-2] == base2[:-2]:
                paireds.append(True)
                del fastqs[idx+1]
            else:
                paireds.append(False)
            idx += 1
        paireds.append(False)

        return fastqs, paireds


    def has_more_fastq(self):
        return self.fastqi + 1 < len(self.fastqs)

    def next_fastq(self):
        n = len(self.fastqs)

        while True:
            self.fastqi += 1
            if self.fastqi >= n:
                return False

            try:
                self.fastq = Fastq(self.fastqs[self.fastqi],
                        paired=self.paireds[self.fastqi])
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

        if not self.selected_testsuites:

            selection = self.selector.selection()
            if not selection:
                tkMessageBox.showerror('no testsuite selected',
                        'please select at least one testsuite before scanning')
                if self.selector.closed:
                    self.selector = TestsuiteSelector(self.testsuite_paths)
                return

            difference = set(selection) - set(self.testsuites.keys())
            self.selector.disable()
            testsuites = load_testsuites(self.testsuite_paths, difference)
            self.testsuites.update(testsuites)
            for name in selection:
                self.selected_testsuites[name] = self.testsuites[name]

        if not self.running:
            # "start"
            self.analyser = analyse.Analyser()

            
            engine.config(**config_params(self.settings.config, self.fastq))

            self.at = AnalyseThread(self.analyser, self.fastq,
                    self.selected_testsuites)
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
            #self.at.join() # deadlocks...!?
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
                lo.error('could not scan %s : %s' %
                        (self.fastq.fname, str(self.at.exception)))
                tkMessageBox.showerror('could not scan',
                        'the scanning of the file "%s" could not be completed : %s' % (
                        self.fastq.fname, self.at.exception))
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
        explorer = JsonExplorer(self.analyser,
                testsuites=self.testsuites, testsuite_paths=self.testsuite_paths)

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
        #self.selector.destroy()
        if self.running:
            self.at.stop()
            self.at.join()
            lo.removeHandler(self.log_handler)


if __name__=='__main__':
            
    from kvarq.gui.settings import Settings
    win = SimpleGUI(Settings())
    tk.mainloop()

