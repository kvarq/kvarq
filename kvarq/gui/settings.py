
import Tkinter as tk
import tkFont
import tkFileDialog
import tkMessageBox
import glob
import os.path

from kvarq.gui.util import open_help, ThemedTk, askopenfilename, get_root_path
from kvarq.genes import load_testsuite, TestsuiteLoadingException
from kvarq.log import lo


class Settings:

    def __init__(self, config):
        self.config = config

    def show(self):
        self.win = ThemedTk(title='settings')

        self.win.rowconfigure(1, weight=1)
        self.win.columnconfigure(1, weight=1)

        row = 0

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
        self.win.bind('<Escape>', lambda x: self.win.close())

    def launch_help(self, e=None):
        open_help(page='gui', anchor='configuring-kvarq')

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

