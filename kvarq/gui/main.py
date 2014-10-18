
from kvarq import VERSION
from kvarq.log import lo, logfn
from kvarq.gui.settings import Settings
from kvarq.gui.simple import SimpleGUI
from kvarq.gui.explorer import DirectoryExplorer
from kvarq.gui.util import open_help, ThemedTk
from kvarq.testsuites import discover_testsuites, load_testsuites
from kvarq.config import default_config

import sys
import Tkinter as tk
import tkMessageBox
import tkFont
import webbrowser
import logging


class GuiLogHandler(logging.Handler):

    def __init__(self, text, scrollfn):
        logging.Handler.__init__(self)
        self.setLevel(logging.DEBUG)

        boldfont = tkFont.Font(text, family='Courier New', size=13, weight='bold')

        text.tag_config('debug', foreground='#888')
        text.tag_config('info', foreground='#080')
        text.tag_config('bold', font=boldfont)
        text.tag_config('warning', background='red', foreground='white')
        text.tag_config('error', background='red', foreground='white')
        text.configure(state='disabled')

        self.text = text
        self.scrollfn = scrollfn
        self.fmt1 = logging.Formatter('[%(asctime)s] -%(levelname)s- %(filename)s:%(lineno)d(%(funcName)s) :: %(message)s')
        self.fmt2 = logging.Formatter('[%(levelname)s] %(message)s')

    def emit(self, record):

        try:
            self.text.insert('end', '')
        except:
            # window destroyed...
            return

        msg = self.fmt2.format(record)
        self.text.config(state='normal')

        if msg[:7] == '[DEBUG]':
            self.text.insert('end', msg + '\n', ('debug',))
        elif msg[:6] == '[INFO]':
            self.text.insert('end', msg[:6], ('info',))
            self.text.insert('end', msg[6:] + '\n')
        elif msg[:9] == '[WARNING]':
            self.text.insert('end', msg[:9], ('warning',))
            self.text.insert('end', msg[9:] + '\n', ('bold',))
        elif msg[:7] == '[ERROR]':
            self.text.insert('end', msg[:7], ('error',))
            self.text.insert('end', msg[7:] + '\n', ('bold',))
        else:
            self.text.insert('end', msg + '\n')

        self.text.configure(state='disabled')
        self.scrollfn()


class MainGUI(ThemedTk):

    def __init__(self, testsuite_paths):
        ThemedTk.__init__(self)

        self.settings = Settings(default_config)
        self.testsuite_paths = testsuite_paths
        self.testsuites = {}

        frame = tk.Frame(self)

        self.scan = tk.Button(frame, text='scan .fastq files', command=self.do_scan)
        self.scan.pack()

        self.explore = tk.Button(frame, text='explore .json files', command=self.do_explore)
        self.explore.pack()

        dummy = tk.Label(frame)
        dummy.pack()

        self.config = tk.Button(frame, text='settings', command=self.do_config)
        self.config.pack()

        self.help = tk.Button(frame, text='help', command=open_help)
        self.help.pack()

        if logfn:
            self.showlog = tk.Button(frame, text='show log file', command=self.do_showlog)
            self.showlog.pack()

        frame.pack(side='left', padx=10)

        outer = tk.Frame(self, borderwidth=1, relief='ridge')
        outer.pack(side='left', expand=True, fill='both', padx=5, pady=5)
        label = tk.Label(outer, text='kvarq log output')
        label.pack()
        frame = tk.Frame(outer)
        frame.pack(expand=True, fill='both')
        self.text = tk.Text(frame) #, state=tk.DISABLED)
        self.text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll = tk.Scrollbar(frame, command=self.text.yview)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.text.config(yscrollcommand=scroll.set)
        self.text.yscrollbar = scroll

        def sys_exit():
            if tkMessageBox.askyesno('quit KvarQ', 'really want to exit KvarQ and close all windows?'):
                sys.exit(0)
        self.protocol("WM_DELETE_WINDOW", sys_exit)

        self.log_handler = GuiLogHandler(self.text, self.scrolldown)
        lo.addHandler(self.log_handler)
        lo.debug('GUI started')

    def do_selection(self, e=None):
        pass

    def do_config(self, e=None):
        self.settings.show()

    def do_scan(self, e=None):
        SimpleGUI(self.settings,
                testsuites=self.testsuites, testsuite_paths=self.testsuite_paths)

    def do_explore(self, e=None):
        DirectoryExplorer(None,
                testsuites=self.testsuites, testsuite_paths=self.testsuite_paths)

    def do_showlog(self, e=None):
        logwin = ThemedTk(title='contents of logfile (%s)' % logfn, geometry=(-200, -200))

        frame = tk.Frame(logwin)
        text = tk.Text(frame)
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll = tk.Scrollbar(frame, command=text.yview)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        text.config(yscrollcommand=scroll.set)
        text.yscrollbar = scroll
        frame.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        for line in file(logfn).readlines():
            text.insert(tk.END, line)
        yy = text.yscrollbar.get()
        text.yscrollbar.set(1-yy[1]+yy[0], 1)
        text.yview('moveto', 1-yy[1]+yy[0])

    def scrolldown(self):
        # autoscrolldown
        yy = self.text.yscrollbar.get()
        self.text.yscrollbar.set(1-yy[1]+yy[0], 1)
        self.text.yview('moveto', 1-yy[1]+yy[0])


if __name__=='__main__':

    win = MainGUI(testsuite_paths=discover_testsuites())
    tk.mainloop()

