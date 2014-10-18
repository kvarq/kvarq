
from kvarq import DOWNLOAD_URL
from kvarq.log import lo, set_debug
from kvarq.analyse import Analyser, DecodingException
from kvarq.analyse import VersionConflictException, TestsuiteVersionConflictException
from kvarq.util import get_root_path, JsonSummary
from kvarq.gui.util import open_help, ThemedTk, BackgroundJob, askopenfilename
from kvarq.gui.tkplot import CoverageWindow, ReadlengthWindow, HitHistogramWindow, \
        MeanCoverageWindow, SpoligoWindow
from kvarq.genes import SNP, TemplateFromGenome
from kvarq.testsuites import update_testsuites

import Tkinter as tk, tkMessageBox, tkFileDialog
import tkFont
import os.path
import glob
import json, codecs
import re
import sys


class DirectoryExplorer:

    def __init__(self, dname, testsuites, testsuite_paths):

        self.testsuites = testsuites
        self.testsuite_paths = testsuite_paths

        if dname:
            self.dname = os.path.abspath(dname)
            self.jpaths = glob.glob(os.path.join(self.dname, '*.json'))
        else:
            jpaths = askopenfilename(
                    initialdir=os.getcwd(),
                    title='Choose .json files to explore',
                    multiple=True,
                    filetypes=[('json files', '*.json')])

            if not jpaths:
                return
            if len(jpaths) == 1:
                self.do_open_json(jpaths[0])
                return

            self.jpaths = [os.path.abspath(jpath) for jpath in jpaths]
            self.dname = os.path.dirname(self.jpaths[0])

        self.win = ThemedTk(title='explore .json files', esc_closes=True,
                geometry=(400, 800))

#        self.menu = tk.Menu(self.win)
#        filemenu = tk.Menu(self.menu)
#        self.menu.add_cascade(label='KvarQ', menu=filemenu)
#        filemenu.add_command(label='Help', command=open_help)
#        self.win.config(menu=self.menu)
#
        self.win.columnconfigure(0, weight=1)
        self.win.columnconfigure(1, weight=0)
        self.win.rowconfigure(0, weight=0)
        self.win.rowconfigure(1, weight=1)
        self.win.rowconfigure(2, weight=0)

        self.dlabel = tk.Label(self.win, text='(no directory chosen)')
        self.dlabel.grid(row=0, column=0, columnspan=2, sticky='ew')

        self.yscroll = tk.Scrollbar(self.win, orient=tk.VERTICAL)
        self.yscroll.grid(row=1, column=1, sticky='ns')
        self.jlist = tk.Listbox(self.win, yscrollcommand=self.yscroll.set)
        self.jlist.grid(row=1, column=0, sticky='nsew')
        self.yscroll["command"] = self.jlist.yview
        self.jlist.bind("<Double-Button-1>", self.open_json)
        self.jlist.bind("<Return>", self.open_json)

        # convert .xls button
        self.convert = tk.Button(self.win, text='summarize...', command=self.summarize)
        self.convert.grid(row=2, column=0, sticky='ew')

        self.update()

        self.jlist.activate(0)
        self.jlist.selection_set(0)
        self.jlist.focus_set()

    def update(self):
        if len(self.dname)>30:
            self.dlabel.config(text='...'+self.dname[-27:])
        else:
            self.dlabel.config(text=self.dname)

        self.jlist.delete(0, tk.END)
        for jpath in self.jpaths:
            self.jlist.insert(tk.END, os.path.basename(jpath))

    def open_json(self, x=None):
        idxs = self.jlist.curselection()
        if not idxs:
            lo.warning('cannot open JsonExplorer : idxs='+str(idxs))
            return
        jpath = self.jpaths[int(idxs[0])]
        self.do_open_json(jpath)

    def do_open_json(self, jpath):
        try:
            je = JsonExplorer(jpath, self.testsuites, self.testsuite_paths)
        except DecodingException, e:
            more = ''
            if isinstance(e, TestsuiteVersionConflictException):
                more += '\n\nYou must load compatible versions of testsuites when ' \
                        'exploring a .json file; KvarQ ships with some testsuite ' \
                        'legacy versions (in the respectevite "legacy/" directory). ' \
                        'Or you may try to find old versions online at ' + DOWNLOAD_URL
            if isinstance(e, VersionConflictException):
                more += '\n\nSome old versions of KvarQ used a file format that ' \
                        'cannot be parsed with this version; you can download older ' \
                        'versions of KvarQ at ' + DOWNLOAD_URL
            tkMessageBox.showerror(
                    'file format error',
                    'cannot load file %s : %s'%(jpath, e)
                    + more)
#        except Exception, e:
#            tkMessageBox.showerror(
#                    'unexpected error',
#                    'cannot load file %s : %s'%(jpath, e))

    def summarize(self, x=None):

        fname = os.path.join(self.dname,  'results.csv')
        i = 2
        while os.path.exists(fname):
            fname = os.path.join(self.dname,  'results%d.csv' % i)
            i += 1

        bj = BackgroundJob('exporting data...')

        self.convert.config(state='disabled')
        text = self.convert.config('text')[4]
        stats = dict(n=0)

        def do_export():

            js = JsonSummary()

            for jpath in self.jpaths:

                if bj.canceled:
                    break

                bj.message = 'extracting from ' + os.path.basename(jpath)

                try:
                    js.add(jpath)
                    stats['n'] += 1
                except Exception as e:
                    lo.error('could not load %s : %s'%(jpath, e))
                    continue

            try:
                js.dump(file(fname, 'w'))
            except IOError as e:
                lo.error('could not write to file %s : %s' % (fname, e))

        def export_done():

            self.convert.config(state='normal')
            self.convert.config(text=text)

            tkMessageBox.showinfo(title='created .xls',
                    message='successfully extracted informations from %d .jsons '
                            'and saved to %s'%(stats['n'], fname))

        bj.start(do_export, export_done)


class JsonExplorer:

    # whishlist
    #   - file info
    #     - version, size, scantime
    #     - readlengths_hist=f(Amin), fastq_type
    #   - lineage
    #     - scores
    #   - resistances
    #     - analyse : (non) synonymous

    def __init__(self, jpath_or_analyser, testsuites, testsuite_paths):
        self.win = ThemedTk(title='json explorer', esc_closes=True,
                geometry=(-200, -200))

        if sys.platform == 'win32':
            self.win.wm_iconbitmap(bitmap = get_root_path('res', 'TPH_DNA.ico'))

        self.win.columnconfigure(0, weight=1)
        self.win.columnconfigure(1, weight=0)
        self.win.rowconfigure(0, weight=0)
        self.win.rowconfigure(1, weight=0)
        self.win.rowconfigure(2, weight=4)

        if isinstance(jpath_or_analyser, Analyser):
            self.analyser = jpath_or_analyser
            name = os.path.basename(self.analyser.fastq.fname)
        else:
            try:

                data = json.load(file(jpath_or_analyser))
                update_testsuites(testsuites, data['info']['testsuites'], testsuite_paths)

                self.analyser = Analyser()
                self.analyser.decode(testsuites, data)
                self.analyser.update_testsuites()
            except Exception, e:
                exc_info = sys.exc_info()
                self.win.destroy()
                raise exc_info[1], None, exc_info[2]
            name = os.path.basename(jpath_or_analyser)

        self.dlabel = tk.Label(self.win, text=name)
        self.dlabel.grid(row=0, column=0, columnspan=2, sticky='ew')

        self.menu = tk.Menu(self.win)
        filemenu = tk.Menu(self.menu)
        self.menu.add_cascade(label='KvarQ', menu=filemenu)
        filemenu.add_command(label='Help', command=open_help)
        self.win.config(menu=self.menu)

        # list of analyses
        self.yscroll1 = tk.Scrollbar(self.win, orient=tk.VERTICAL)
        self.yscroll1.grid(row=1, column=1, sticky='ns')
        self.alist = tk.Listbox(self.win, height=len(self.analyser.testsuites)+1, yscrollcommand=self.yscroll1.set)
        self.alist.grid(row=1, column=0, sticky='nsew')
        self.yscroll1["command"]  =  self.alist.yview
        self.alist.bind("<Double-Button-1>", self.show_analyses)
        self.alist.bind("<Return>", self.show_analyses)

        # fill in list of analyses, prepare dict of coverages
        self.anames = ['info']
        self.alist.insert(tk.END, 'info')
        for aname, testsuite in self.analyser.testsuites.items():
            self.anames.append(aname)
            result = self.analyser.results[aname]
            if type(result)==list:
                result = '; '.join(result)
            self.alist.insert(tk.END, aname + ': ' + result)

        # list of coverages
        self.yscroll2 = tk.Scrollbar(self.win, orient=tk.VERTICAL)
        self.yscroll2.grid(row=2, column=1, sticky='ns')
        self.clist = tk.Listbox(self.win, yscrollcommand=self.yscroll2.set)
        self.clist.grid(row=2, column=0, sticky='nsew')
        self.yscroll2["command"]  =  self.clist.yview
        self.clist.bind("<Double-Button-1>", self.show_coverage)
        self.clist.bind("<Return>", self.show_coverage)

        self.current = None
        self.alist.activate(0)
        self.alist.selection_set(0)
        self.alist.focus_set()

        self.after_id = None
        def close_win(a=None):
            if self.after_id:
                self.win.after_cancel(self.after_id)
            self.win.destroy()
        self.win.close = close_win
        self.win.protocol('WM_DELETE_WINDOW', close_win)
        self.poll()

    def show_analyses(self, x):

        idxs = self.alist.curselection()
        if not idxs: return
        aname = self.anames[int(idxs[0])]
        if aname == 'spoligo':
            spoligos = [-1] * 43
            for test in self.analyser.testsuites['spoligo'].tests:
                spoligos[test.genotype.number] = self.analyser[test].mean()
            SpoligoWindow(spoligos)
#            fw = FigureWin()
#            self.illustrator.plot_spoligo(fw.plt)
#            fw.win.focus_force()

    def show_coverage(self, x):

        idxs = self.clist.curselection()
        if not idxs: return

        # special info 'analysis'
        if self.aname == 'info':
            iname = self.infos[int(idxs[0])]
            if iname == 'readlengths...':
                ReadlengthWindow(self.analyser.stats['readlengths'])
#                fw = FigureWin()
#                for rl, n in enumerate(self.analyser.stats['readlengths']):
#                    fw.plt.plot([rl, rl], [0, n], 'k')
#                fw.win.focus_force()
            if iname == 'mean coverage...':
#                coverages = [
#                        n/len(self.analyser[i])
#                        for i, n in enumerate(self.analyser.stats['nseqbasehits'])]
                mean_coverages = [coverage.mean(include_margins=False)
                        for coverage in self.analyser.coverages.values()]
                MeanCoverageWindow(mean_coverages)
#                fw = FigureWin()
#                fw.plt.hist(coverages)
#                fw.win.focus_force()
            if iname == 'hits/template...':
                nseqhits = self.analyser.stats['nseqhits']
                HitHistogramWindow([
                    sum([nseqhits[idx] for idx in self.analyser.get_indexes(coveragename)])
                    for coveragename in self.analyser.coverages])
#                fw = FigureWin()
#                for rl, n in enumerate(self.analyser.stats['nseqhits']):
#                    fw.plt.plot([rl, rl], [0, n], 'k')
#                fw.win.focus_force()
            return

        test = self.tests_sorted[int(idxs[0])]
        try:
            coverage = self.analyser[test]
        except KeyError, e:
            tkMessageBox.showinfo(
                    title='test not found',
                    message='"%s" not found in .json' % str(test))
            return
        cw = CoverageWindow(test, coverage)
#        fw = FigureWin()
#        self.illustrator.plot_coverage(fw.plt, self.aname, test)
#        fw.win.focus_force()

    def update(self):
        idxs = self.alist.curselection()
        if not idxs: return
        self.aname = self.anames[int(idxs[0])]
        self.clist.delete(0, tk.END)

        # special info 'analysis'
        if self.aname == 'info':
            self.infos = [
                    'fastq : ' + ', '.join(self.analyser.fastq_filenames),
                    'size : ' + ', '.join(['%.2f MB'%(fastq_size/1024.**2)
                            for fastq_size in self.analyser.fastq_sizes]),
                    'readlength : %d'%self.analyser.fastq_readlength,
                    'records_approx : %s'%str(self.analyser.fastq_records_approx or '?'),
                    'scantime : %d s'%int(self.analyser.scantime),
                    'config : ' + ' '.join([str(k)+'='+str(v)
                            for k, v in self.analyser.config.items()]),
                    '',
                    ## stats
                    'readlengths...',
                    'mean coverage...',
                    'hits/template...',
                    'records_parsed : %d'%self.analyser.stats.get('records_parsed', -1),
                    'progress : %.1f %%'%(float(self.analyser.stats['progress'])*100),
                ]

            for info in self.infos:
                self.clist.insert(tk.END, info)
            return

        # src : http://stackoverflow.com/questions/5254021/
        #       python-human-sort-of-numbers-with-alpha-numeric-but-in-pyqt-and-a-lt-oper
        def _human_key(key):
            parts = re.split('(\d*\.\d+|\d+)', str(key))
            return tuple((e.swapcase() if i % 2 == 0 else float(e))
                    for i, e in enumerate(parts))

        tests = self.analyser.testsuites[self.aname].tests
        self.tests_sorted = sorted(tests, key=_human_key)

        for i, test in enumerate(self.tests_sorted):
            try:
                coverage = self.analyser[test]
            except KeyError, e:
                self.clist.insert(tk.END, '(test %s not found in .json)' % str(test))
                continue
            seqmean = coverage.seqmean()
            mean = coverage.mean(include_margins=False)

            sign = ''
            if coverage.mixed():
                sign += '~'
            if isinstance(test.template, TemplateFromGenome) and \
                    not isinstance(test.template, SNP):
                sign += '+' * len(test.template.mutations(coverage))
            else:
                if test.template.validate(coverage):
                    sign += '+'

            hits = ''
            if 'nseqhits' in self.analyser.stats:
                idxs = self.analyser.get_indexes(test)
                nseqhits = self.analyser.stats['nseqhits']
                hits = '%d hits ' % sum([nseqhits[idx] for idx in idxs])

            self.clist.insert(tk.END, '%s %s %s(mean %.1f/%.1f)'%(
                    sign, str(test), hits, seqmean, mean))

    def poll(self):
        now = self.alist.curselection()
        if now != self.current:
            self.update()
            self.current = now
        self.after_id = self.win.after(250, self.poll)



if __name__ == '__main__':

    import sys, argparse

    parser = argparse.ArgumentParser(description='''
            interactively explore .json files in directory previously
            generated by "kvarq scan" ''')

    parser.add_argument('-d', '--debug', action='store_true',
            help='output log information at a debug level')
    parser.add_argument('directory', nargs='?',
            help='directory to explore')

    args = parser.parse_args(sys.argv[1:])

    if args.debug:
        set_debug()

    win = DirectoryExplorer(args.directory, testsuites={})
    tk.mainloop()

