
import Tkinter as tk
import tkFont
import tkFileDialog
import tkMessageBox
import glob
import os.path

from kvarq.gui.util import open_help, ThemedTk, askopenfilename, get_root_path
from kvarq.genes import load_testsuite, TestsuiteLoadingException
from kvarq.config import default_config
from kvarq.log import lo


class Settings:

    def __init__(self, config=default_config, selection=None, testsuites=None):
        self.config = config
        if testsuites is None:
            testsuites = dict()
        self.testsuites = testsuites
        if selection is None:
            selection = testsuites.keys()
        self.selection = selection

    def get_testsuites(self):
        ''' returns a testsuites dictionary based on current selection '''
        return dict([(name, self.testsuites[name]) for name in self.selection])

    def show(self):
        self.win = ThemedTk(title='settings')

        self.win.rowconfigure(1, weight=1)
        self.win.columnconfigure(1, weight=1)

        row = 0

        l = tk.Label(self.win, text='Testsuites', font=self.win.boldfont)
        l.grid(row=row, column=0, sticky='nw')

        self.testsuites_list = tk.Listbox(self.win, height=8, selectmode='multiple')
        self.testsuites_list.grid(row=row, rowspan=2, column=1, columnspan=2, sticky='nesw')
        for name in self.testsuites.keys():
            self.testsuites_list.insert('end', name)
            if name in self.selection:
                self.testsuites_list.selection_set('end')

        row += 1

        l = tk.Label(self.win, text='select from list')
        l.grid(row=row, column=0, sticky='nw')

        row += 1

        self.testsuites_load = tk.Button(self.win, text='load...',
                command=self.load_testsuites)
        self.testsuites_load.grid(row=row, column=2)

        row += 1
        self.win.rowconfigure(row, minsize=10)
        row += 1

        l = tk.Label(self.win, text='Engine configuration', font=self.win.boldfont)
        l.grid(row=row, column=0, columnspan=2, sticky='w')

        self.show_help = tk.Button(self.win, text='?',
                command=self.launch_help)
        self.show_help.grid(row=row, column=2, sticky='e')

        self.entries = {}
        for name in self.config.keys():
            row += 1
            l = tk.Label(self.win, text=name)
            l.grid(row=row, column=0, sticky='w')
            self.entries[name] = tk.Entry(self.win)
            self.entries[name].grid(row=row, column=1, columnspan=2, sticky='ew')
            self.entries[name].insert(0, self.config[name])

        row += 1
        self.win.rowconfigure(row, minsize=10)
        row += 1
        frame = tk.Frame(self.win)
        frame.grid(row=row, column=0, columnspan=3)
        save = tk.Button(frame, text='save', command=self.save_cb)
        save.pack(side='left')
        cancel = tk.Button(frame, text='cancel', command=self.cancel_cb)
        cancel.pack(side='left')

    def launch_help(self, e=None):
        open_help(page='gui.html', anchor='configuring-kvarq')

    def load_testsuites(self, e=None):

        initialdir = get_root_path('testsuites')
        if not os.path.exists(initialdir):
            initialdir = get_root_path()

        fnames = askopenfilename(
                title='open additional testsuites',
                initialdir=initialdir,
                multiple=True,
                filetypes=[('testsuites', '*.py')])

        replaced = []

        for fname in fnames:
            try:
                name, testsuite = load_testsuite(fname)
                if name in self.testsuites:
                    replaced.append('testsuite "%s" %s -> %s'%(
                            (name, self.testsuites[name].version, testsuite.version)))
                    lo.info('replaced testsuite "%s" %s -> %s'%(
                            (name, self.testsuites[name].version, testsuite.version)))
                else:
                    self.testsuites_list.insert('end', name)
                    self.testsuites_list.selection_set('end')
                    self.selection.append(name)
                    lo.info('loaded testsuite "%s" %s' % (name, testsuite.version))
                self.testsuites[name] = testsuite
            except TestsuiteLoadingException, e:
                tkMessageBox.showerror(
                        'cannot load testsuite',
                        'testsuite "%s" cannot be loaded : %s' % (fname, e))

        if replaced:
            tkMessageBox.showinfo('replacements', '\n'.join(replaced))

    def save_cb(self, e=None):

        # check whether parameters have correct format
        newconfig = dict(self.config)
        for name in self.config.keys():
            try:
                newconfig[name] = int(self.entries[name].get())
            except ValueError, e:
                tkMessageBox.showerror(
                        'invalid value',
                        'parameter "%s" must be an integer' % name)
                return

        # save new settings
        for name, value in newconfig.items():
            self.config[name] = value

        while self.selection: self.selection.pop()
        for i in range(self.testsuites_list.size()):
            if self.testsuites_list.selection_includes(i):
                self.selection.append(self.testsuites_list.get(i))

        self.win.destroy()

    def cancel_cb(self, e=None):
        self.win.destroy()



if __name__ == '__main__':

    settings = Settings()
    settings.show()
    tk.mainloop()
    # settings.config changes depending what button was clicked
    settings.show()
    tk.mainloop()

