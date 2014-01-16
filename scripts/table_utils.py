
from kvarq.util import csv_xls_writer


class csv_xls_wrapper:

    class xls_iter_wrapper:
        def __init__(self, ws):
            self.ws = ws
            self.row = -1
        def __iter__(self):
            return self
        def next(self):
            self.row += 1
            if self.row == self.ws.nrows:
                raise StopIteration
            return [
                    self.ws.cell(self.row, col).value
                    for col in range(self.ws.ncols)
                ]

    def __init__(self, fname):
        if fname.endswith('.csv'):
            import csv
            self.iter = csv.reader(file(fname)).__iter__()
        elif fname.endswith('.xls'):
            try:
                import xlrd
            except ImportError:
                raise ImportError('could not locate package "xlrd" needed to parse .xls file')
            ws = xlrd.open_workbook(fname).sheet_by_index(0)
            self.iter = csv_xls_wrapper.xls_iter_wrapper(ws)
        else:
            raise Exception('can only read .csv and .xls files')

    def __iter__(self):
        return self.iter


