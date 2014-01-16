"""
GUI clients
"""

# src : http://stackoverflow.com/questions/8691655/python-put-window-on-top-tkinter-pyobjc
def lift_window(win):
    win.lift()
    win.call('wm', 'attributes', '.', '-topmost', True)
    def after_lift():
        win.call('wm', 'attributes', '.', '-topmost', False)
    win.after_idle(after_lift)

