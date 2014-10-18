
import unittest

from kvarq.util import TextHist


class UtilTest(unittest.TestCase):

    def test_texthist(self):
        # min bin width is 1
        data = [0, 0, 1]
        hist = TextHist().draw(sorted(data), indexed=False)
        assert '100%' in hist
        data = [0, 0, 2]
        hist = TextHist().draw(sorted(data), indexed=False)
        assert '66%' in hist
        # all zero
        data = [0, 0, 0]
        hist = TextHist().draw(sorted(data), indexed=False)
        print(hist)
        assert 'CANNOT' in hist
        # zero width bin
        data = [1, 1, 1]
        hist = TextHist().draw(sorted(data), indexed=False)
        assert 'CANNOT' in hist
        # no data
        data = []
        hist = TextHist().draw(sorted(data), indexed=False)
        assert 'CANNOT' in hist


if __name__ == '__main__': unittest.main()

