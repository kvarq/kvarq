
from kvarq import VERSION, DOC_URL
from kvarq.log import lo
from kvarq.util import is_exe, is_app, get_root_path, get_help_path

import os, sys, os.path
import webbrowser
import threading

import Tkinter as tk
import tkFont
import tkFileDialog

class ThemedTk(tk.Tk):

    def __init__(self, title=None, esc_closes=False, geometry=None):
        ''' :param esc_closes: whether hitting the ``<Escape>`` key should close
                the window (via overridable :py:meth:`close` callback)

            :param geometry: tuple ``(width, height)`` where both parameters can
                be either pixels or fractions (value between 0-1) or difference
                from total screen dimension (value below 0) '''

        tk.Tk.__init__(self)
        if title:
            title = ' -- ' + title
        else:
            title = ''
        self.wm_title('KvarQ ' + str(VERSION) + title)

        if sys.platform == 'win32':
            self.wm_iconbitmap(bitmap = get_root_path('res', 'TPH_DNA.ico'))

        self.monospace = tkFont.Font(self, family='Courier New', weight=tkFont.BOLD)
        self.boldfont = tkFont.Font(self, weight='bold')

        if esc_closes:
            self.bind('<Escape>', lambda x: self.close())

        if geometry:
            sw = self.winfo_screenwidth()
            sh = self.winfo_screenheight()

            w, h = geometry
            if w<0: w = sw + w
            elif w<1: w = int(sw * w)
            if h<0: h = sh + h
            elif h<1: h = int(sh * h)

            w = max(200, min(w, sw - 100))
            h = max(200, min(h, sh - 200))
            self.geometry('%dx%d+%d+%d'%(w,h,50,(sh-h)/2))

        self.focus_force()

    def close(self, x=None):
        self.destroy()


def focus_force():
    if sys.platform == 'darwin':
        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')


def open_help(page='index', anchor=None):
    # unfortunately, anchors don't work under windows:
    # http://stackoverflow.com/questions/6374761/python-webrowser-open-url-with-bookmarks-like-www-something-com-file-htmltop
    webbrowser.open(get_help_path(page, anchor, need_url=True))


class BackgroundJob(tk.Tk):

    def __init__(self, title):
        tk.Tk.__init__(self)
        self.title(title)
        self.label = tk.Label(self)
        self.label.pack(expand=1, fill='x')
        self.cancel = tk.Button(self, text='cancel', command=self.cancel_cb)
        self.cancel.pack()
        self.geometry('300x150')
        self.resizable(0, 0)

        # will be set if "cancel" is clicked
        self.canceled = False
        # will be used in GUI thread to set label value
        self.message = ''
        # use this for data transfer between threads
        self.data = None

    def start(self, run, done_cb=None):
        self.thread = threading.Thread()
        self.thread.run = run
        self.thread.start()
        self.done_cb = done_cb
        self.update()

    def update(self):
        if self.thread.is_alive():
            self.label.config(text=self.message)
            self.after(100, self.update)
        else:
            self.thread = None
            self.destroy()
            if self.done_cb:
                self.done_cb()

    def cancel_cb(self, e=None):
        self.canceled = True

# src : http://stackoverflow.com/questions/4116249/parsing-the-results-of-askopenfilenames
def askopenfilename(*a, **b):
    ret = tkFileDialog.askopenfilename(*a, **b)
    if type(ret) == tuple:
        return ret
    tmp = tk.Tk()
    tmp.withdraw()
    return tmp.tk.splitlist(ret)

