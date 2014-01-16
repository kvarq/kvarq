'''
High level interface to the scanning process.  The :py:class:`Analyser` is
the central instance to communicate with the :py:mod:`engine` and load/save the
results to/from ``.json`` files.
'''

from kvarq import VERSION
from kvarq.log import lo, tictoc, format_traceback
from kvarq import engine
from kvarq import genes
from kvarq.util import TextHist
from kvarq.fastq import Fastq
from kvarq.legacy import convert_legacy_data

import json, codecs
import time
import os.path
from distutils.version import StrictVersion
from collections import Counter, OrderedDict
import re
import sys


class Coverage:
    '''
    This class applies :py:class:`kvarq.engine.Hit` to a
    :py:class:`kvarq.genes.Sequence`, keeping track of matching and non
    matching bases.
        
    A coverage does not record the position inside the genome (this is stored
    in the corresponding :py:class:`kvarq.genes.Template` class) but it keeps
    track of the margins on both sides that were used during scanning.
        
    example: if SNP1000 is scanned with a margin=10 then the coverage will
    have ``len(coverage.coverage)=21``, ``coverage.start=10``,
    ``coverage.stop=11`` with ``coverage.coverage[10]`` corresponding to base
    pair 1000 on the genome
    '''

    def __init__(self, plus_seq):
        '''
        :param plus_seq: a :py:class:`kvarq.genes.Sequence` on the ``+`` strand
            to which hits will be mapped
        '''
        self.plus_seq = plus_seq
        self.minus_seq = plus_seq.reverse()
        self.coverage = [0] * len(plus_seq)
        self.mutations = {}
        #TODO rename to left, right
        self.start = plus_seq.left
        self.stop = len(plus_seq) - plus_seq.right

    def apply_hit(self, hit, hitseq, on_plus_strand):
        '''
        :param hit: a :py:class:`kvarq.engine.Hit` to be mapped on the sequence
        :param hitseq: string of bases found in the ``.fastq`` file
        :param on_plus_strand: whether the hit refers to the ``+`` strand

        applies a hit to the sequence
        '''
        if on_plus_strand:
            seq = self.plus_seq
        else:
            seq = self.minus_seq

        start = max(0, hit.seq_pos)
        stop = start+hit.length

        for i, j in enumerate(range(start, stop)):
            c_j = seq.plus_idx(j)
            c_b = seq.plus_base(hitseq[i])
            self.coverage[c_j] += 1
            if hitseq[i] != seq[j]:
                self.mutations[c_j] = self.mutations.get(c_j, '') + c_b

    def bases_at(self, idx):
        ''' :returns: array of ``bases=>coverage`` at specified position '''
        m = self.mutations.get(idx, '')
        ret = { self.plus_seq[idx]:self.coverage[idx]-len(m) }
        for b in set(m):
            ret[b] = m.count(b)
        return ret

    def mean(self, include_margins=True):
        '''
        :param include_margins: whether to include the two flanks
        :returns: the average depth of coverage (optionally including the flanks)
        '''
        if include_margins:
            return sum(self.coverage)/float(len(self.coverage))
        else:
            return sum(self.coverage[self.start:self.stop])/float(self.stop-self.start)

    def std(self, include_margins=True):
        '''
        :param include_margins: whether to include the two flanks
        :returns: the average depth of coverage (optionally including the flanks)
        '''
        m = self.mean()
        if include_margins:
            xs = self.coverage
        else:
            xs = self.coverage[self.start:self.stop+1]
        return (sum([(x-m)**2 for x in self.coverage])/len(self.coverage))**.5

    def seqmean(self):
        ''' :returns: mean coverage of sequence, *not* counting mutations '''
        seq = [self.coverage[i] - len(self.mutations.get(i, ''))
                for i in range(self.start, self.stop)]
        return sum(seq)/float(len(seq))

    def __str__(self):
        ''' :returns: string containing mean, std & serialized coverage '''
        return '(mean %.2f std %.2f) '%(self.mean(), self.std()) \
                + ' ' + self.serialize()

    def serialize(self):
        '''
        :returns: a stringified and human-readable representation of
            ``.coverage`` and ``.mutations`` '''
        cov = '-'.join([str(c) for c in self.coverage])
        mut = '-'.join(['%d[%s]'%(idx, ''.join(sorted(self.mutations[idx])))
                for idx in sorted(self.mutations.keys())])
        return cov + ' ' + mut

    def deserialize(self, serialized_coverage):
        ''' :param serialized_coverage: string as returned from ``.serialize()``

            recreates ``.coverage`` and ``.mutations`` as specified by serialized
            string '''
        c_s, space, m_s = serialized_coverage.partition(' ')
        self.coverage = [int(x) for x in c_s.split('-')]
        if m_s:
            self.mutations = dict([(int(x[:x.index('[')]), x[x.index('[')+1:x.index(']')])
                    for x in m_s.split('-')])
        else:
            self.mutations = {}

    def __len__(self):
        ''' :returns: length of ``.coverage`` '''
        return self.coverage.__len__()

    def __getitem__(self, idx):
        ''' :returns: depth of coverage at specified position '''
        return self.coverage.__getitem__(idx)



class DecodingException(Exception):
    ''' issued when :py:class:`Analyser` cannot be decode()d '''

class VersionConflictException(DecodingException):
    ''' issued when :py:class:`Analyser` cannot be decode()d 
        because the file version is incompatible with the current
        KvarQ version '''

class TestsuiteVersionConflictException(DecodingException):
    ''' issued when :py:class:`Analyser` cannot be decode()d 
        because the testsuites that were used to generate a ``.json``
        file are not compatible with the version of the currently
        loaded testsuites '''

class DataInconcistencyException(DecodingException):
    ''' issued when :py:class:`Analyser` cannot be decode()d due to some
        inconsistency in the decoded data '''

class Analyser:

    '''
    Main class for analyzing raw genome data. the ``.fastq`` file is first
    scanned for hits, then these hits are used to fill intermediary data
    structures (:py:class:`kvarq.analyse.Coverage`) and finally genetic
    information as defined in :py:mod:`kvarq.genes` is extracted.

    Normal scenario for scanning a ``.fastq`` file and saving results
    to a ``.json`` file::

        analyser = Analyser()
        analyser.scan(fastq, testsuites) # calls engine
        analyser.update_testsuites() # creates .results
        json.dump(analyser.encode(), file(fname, 'w'))

    Alternatively, data can be read from ``.json`` file and results can
    be regenerated without scanning the ``.fastq`` file::

        analyser = Analyser()
        analyser.decode(testsuites, json.load(file(fname)))
        analyser.update_testsuites() # recreate .results

    the compatibility of the data structure generated by calling
    :py:meth:`encode` is defined by :py:data:`kvarq.VERSION` (following
    the usual semantics; see http://semvar.org).  The different testsuites
    define their own compatibility sematics (see ``TESTSUITE_VERSION`` and
    ``GENES_COMPATIBILITY`` module globals).
    '''

    def __init__(self):
        # configuration
        self.config = None

        # other information
        self.fastq_fname = None
        self.fastq_size = None
        self.fastq_readlength = None
        self.fastq_records_approx = None
        self.spacing = None

        # direct results from scanning
        self.hits = None 
        self.stats = None
        self.scantime = 0

        # list of coverages will be generated upeon scanning/decoding
        self.coverages = None

        # final results
        self.results = None


    def load_coverages(self, testsuites, spacing=None):
        '''
        :param spacing: how many bases are added on either side to the
            sequence from templates that are read from a reference genome

        :returns: an :py:class:`collections.OrderdDict` of :py:class:`.Coverage`
            indexed by ``str(Test)`` -- each of these coverages was initialized
            with a :py:class:`kvarq.genes.Sequence` with flanks of the
            length specified by ``spacing``
        '''

        if not spacing and not self.spacing:
            spacing = 25 #TODO determine margins depending on self.fastq
        self.spacing = spacing

        coverages = OrderedDict()

        for name, testsuite in testsuites.items():
            for test in testsuite.tests:
                if isinstance(test.template, genes.DynamicTemplate):
                    seq = test.template.seq(spacing=self.spacing)
                else:
                    seq = test.template.seq()

                coverages[str(test.template)] = Coverage(seq)

        return coverages

    def coverage_at(self, i):
        '''
        returns a :py:class:`Coverage` by index. indices are defined when
        sequences are serialized (see :py:meth:`.load_coverages`). indices can
        exted up to ``2*len(coverages)`` (mapping reverse sequences in
        ``Hit.seq_nr``)
        '''
        n = len(self.coverages)
        if i >= n:
            i -= n
        return self.coverages[self.coverages.keys()[i]]

    def get_indexes(self, thing):
        '''
        :param thing: integer, string or :py:class:`kvarq.genes.Test`
        :returns: ``[forward_idx, reverse_idx]`` that can be used to access
            coverages by integer -- see :py:meth:`.__getitem__``
        '''
        if isinstance(thing, genes.Test):
            idx1 = self.coverages.keys().index(str(thing.template))
        else:
            idx1 = self.coverages.keys().index(thing)

        # forward & reverse
        return [idx1, idx1 + len(self.coverages)]

    def __len__(self):
        ''' :returns: number of coverages '''
        return len(self.coverages)

    def __getitem__(self, thing):
        '''
        :param thing: can be :py:class:`kvarq.genes.Test`, its string
            representation or and integer. integers are passed on to
            :py:meth:`.coverage_at`
        :returns: :py:class:`.Coverage`
        '''
        if type(thing) == int:
            return self.coverage_at(thing)
        elif isinstance(thing, genes.Test):
            return self.coverages[str(thing.template)]
        else:
            return self.coverages[str(thing)]

    def scan(self, fastq, testsuites, do_reverse=True):
        '''
        :param fastq: :py:class:`kvarq.fastq.Fastq` file to scan
        :param testsuites: dictionary of instances of
            :py:class:`kvarq.genes.Testsuite`

        initiates a :py:func:`kvarq.engine.findseqs` and fills the attributes
        ``.hits``, ``.stats`` and ``.coverages``

        note that the call to :py:func:`kvarq.engine.findseqs` can raise a
        :py:class:`kvarq.engine.FastqFileFormatException`
        '''

        self.fastq = fastq
        self.fastq_fname = fastq.fname
        self.fastq_size = fastq.size
        self.fastq_readlength = fastq.readlength
        self.fastq_records_approx = fastq.records_approx

        self.testsuites = testsuites
        self.coverages = self.load_coverages(testsuites)

        self.config = engine.get_config()

        seqs = [coverage.plus_seq.bases for coverage in self.coverages.values()]
        if do_reverse:
            seqs += [coverage.minus_seq.bases for coverage in self.coverages.values()]

        # do the scanning
        t0 = time.time()
        ret = engine.findseqs(self.fastq.fname, seqs)
        lo.debug('found %d hits' % len(ret['hits']))
        self.stats = ret['stats']
        self.hits = ret['hits']
        self.scantime = time.time() - t0

        self.update_coverages()


    @tictoc('update_coverages')
    def update_coverages(self):
        ''' applies ``.hits`` to ``.coverages``; this method is called from
            :py:meth:`.scan` but can also be called manually after a call to
            :py:meth:`.decode` if for some reason the ``.coverages`` have to
            be regenerated '''

        assert self.hits is not None, 'cannot update coverages without .hits'
        assert self.fastq is not None, 'cannot update coverages without .fastq'

        for hit in self.hits:

            coverage = self.coverage_at(hit.seq_nr)
            hitseq = self.fastq.readhit(hit)
            coverage.apply_hit(hit, hitseq, hit.seq_nr < len(self.coverages))


    def update_testsuites(self):
        ''' creates ``.results`` by calling :py:meth:`kvarq.genes.Testsuites.analyse`
            on every testsuite -- call this method after :py:meth:`.scan` or
            :py:meth:`.decode` '''
        self.results = {}
        for name, testsuite in self.testsuites.items():
            try:
                self.results[name] = testsuite.analyse(self)
            except Exception, e:
                lo.error('testsuite "%s" : %s [%s]' % (
                        name, e, format_traceback(sys.exc_info())))
                self.results[name] = 'ERROR : ' + str(e)

    @tictoc('encode')
    def encode(self, hits=False):
        ''' returns an object containing the following serializable information

            - ``analyses`` : final scanning results
            - ``info`` : meta-information about the file and scanning parameters
            - ``stats`` : scanning statistics
            - ``coverages`` : intermediate results (a :py:class:`Coverage` for
              every :py:class:`kvarq.genes.Test` in every used
              :py:class:`kvarq.genes.Testsuite`)
            - ``hits`` (optional) : direct results form scanning
        '''

        more ={}
        if hits:
            more['hits'] = self.hits

        return dict(
                analyses=self.results,
                info={
                    'format':'kvarq',
                    'fastq':self.fastq_fname,
                    'size':self.fastq_size,
                    'readlength':self.fastq_readlength,
                    'records_approx':self.fastq_records_approx,
                    'scantime':self.scantime,
                    'when':time.asctime(time.localtime()),
                    'version':VERSION,
                    'config':self.config,
                    'spacing':self.spacing,
                    'testsuites':dict([(name, testsuite.version)
                            for name, testsuite in self.testsuites.items()]),
                },
                stats=self.stats,
                coverages=[(name, coverage.serialize())
                    for name, coverage in self.coverages.items()],
                **more
            )


    @tictoc('decode')
    def decode(self, testsuites, data):
        '''
        :param testsuites: dictionary of :py:class:`kvarq.genes.Testsuite`
        :param data: dictionary as returned by :py:meth:`.encode`

        regenerates attributes as they were after the call to :py:meth:`.scan`
        previous to the call to :py:meth:`.encode` that generated the data
        passed as parameter ``data``

        the attributes ``.fastq`` can only be re-generated if the ``.fastq``
        file can still be found under location specified in the ``.json``
        file

        this method uses :py:func:`kvarq.legacy.convert_legacy_data` to load
        data generated by older version of KvarQ. if the testsuites that were
        used to generate the data are not compatible with the testsuites currently
        loaded, a :py:class:`.VersionConflictException` is risen. tests/coverages
        not used by any of the testsuites specified in the parameter are simply
        ignored. testsuites not found in the ``data`` are also disregarded.
        '''

        data = convert_legacy_data(testsuites, data)

        self.config = data['info']['config']
        self.fastq_fname = data['info']['fastq']
        self.fastq_size = data['info']['size']
        self.fastq_readlength = data['info'].get('readlength', -1)
        self.fastq_records_approx = data['info'].get('records_approx', -1)
        self.stats = data['stats']
        self.scantime = data['info'].get('scantime', -1)

        if 'hits' in data:
            self.hits = [engine.Hit(*hit) for hit in data['hits']]
        else:
            self.hits = None

        fastqname = data['info']['fastq']
        if os.path.isfile(fastqname):
            lo.info('found .fastq file : ' + fastqname)
            self.fastq = Fastq(fastqname)
        else:
            lo.info('cannot load .fastq file : ' + fastqname)
            self.fastq = None

        # prepare testsuites
        self.testsuites = {}
        for name, version in data['info']['testsuites'].items():
            if name in testsuites:

                testsuite = testsuites[name]
                json_v = StrictVersion(version)
                kvarq_v = StrictVersion(testsuite.version)

                if json_v > kvarq_v or json_v.version[0] != kvarq_v.version[0]:
                    raise TestsuiteVersionConflictException(
                            'version conflict testsuite "%s" : '
                            '.json version "%s" not compatible with current version "%s"' %
                            (name, version, testsuite.version))

                lo.debug('loading testsuite %s (%s)'%(name, str(testsuite)))
                self.testsuites[name] = testsuite

            else:
                lo.warning('testsuite "%s" not loaded -> ignoring some results in '
                        '.json file' % name)

        # convert array of (coveragename, serialized_coverage) into OrderedDict
        # (sequences read from templates)
        templates = dict()
        for testsuite in testsuites.values():
            for test in testsuite.tests:
                templates[str(test.template)] = test.template

        self.spacing = data['info']['spacing']
        self.coverages = OrderedDict()
        for name, serialized_coverage in data['coverages']:
            if not name in templates:
                # newer versions of testsuites can discard tests (backwards
                # compatible change with increase in minor version number)
                continue
                #raise DecodingException('template "%s" not found in testsuites %s' %
                #        (name, ', '.join(self.testsuites.keys())))

            template = templates[name]
            if isinstance(template, genes.DynamicTemplate):
                seq = template.seq(spacing=self.spacing)
            else:
                seq = template.seq()

            coverage = Coverage(seq)
            coverage.deserialize(serialized_coverage)
            self.coverages[name] = coverage



class AnalyserJson:

    ''' helper class to handle .json files created from
        :py:class:`Analyser`.encode() '''

    def __init__(self, jpath, minver=None):

        ''' loads data from specified .json file ``jpath`` and does some
            file format checking (including minimal version as specified
            by ``minver``) '''

        try:
            self.data = json.load(codecs.open(jpath, encoding='utf-8'))
        except ValueError, e:
            raise DecodingException, 'not valid .json format : '+str(e)

        if not 'info' in self.data:
            raise DecodingException, 'not valid file format : "info" key missing'
        if not 'format' in self.data['info'] or not self.data['info']['format']=='kvarq':
            raise DecodingException, 'not valid file format : "info"/"format" != "kvarq"'

        if minver:
            minver = StrictVersion(minver)
            dataversion = StrictVersion(self.data['info']['version'])
            if dataversion<minver:
                raise VersionConflictException, '.json format too old : %s < %s'%(dataversion, minver)


    @property
    def analyses(self):
        return self.data['analyses'].items()

