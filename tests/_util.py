
import logging

from kvarq.log import lo

# src : http://stackoverflow.com/questions/4398967/python-unit-testing-automatically-running-the-debugger-when-a-test-fails
import functools, pdb, sys
def debug_on(*exceptions):
    if not exceptions:
        exceptions = (AssertionError, )
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                pdb.post_mortem(sys.exc_info()[2])
        return wrapper
    return decorator 

class NeedleHandler(logging.Handler):

    def __init__(self, needle, needle_level):
        logging.Handler.__init__(self, logging.DEBUG)
        self.needle = needle
        self.needle_level = needle_level
        self.found = False

    def handle(self, record):
        if self.needle_level is None or self.needle_level == record.levelno:
            if self.needle in record.msg:
                self.found = True

def lo_exceptor(needle, level=None, suppress=False):
    ''' returns a function which will assert the specified log
        string/level when called; optionally suppresses other
        handlers until asserted '''

    levels = {}
    if suppress:
        for handler in lo.handlers:
            levels[handler] = handler.level
            handler.level = logging.FATAL

    needlehandler = NeedleHandler(needle, level)
    lo.addHandler(needlehandler)

    def lo_assert():
        msg = 'expected log message "%s"' % needle
        if level is not None:
            msg += ' (level %s)' % logging.getLevelName(level)
        assert needlehandler.found, msg

        lo.removeHandler(needlehandler)
        for handler in lo.handlers:
            handler.level = levels[handler]

    return lo_assert

