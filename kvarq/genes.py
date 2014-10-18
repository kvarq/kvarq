'''
defines classes to represent genetic information and interfaces for
testsuites
'''

import kvarq
from kvarq.util import get_root_path
from kvarq.log import lo, format_traceback

import os.path, glob
import sys
import re
from distutils.version import StrictVersion

'''defines which testsuites can be loaded by this version of KvarQ. whenever
new changes are introduced that break backwards-compatibility, the first number
is increased by one. introduction of new features that retain compatibility
with old testsuites will increase the second number by one'''
COMPATIBILITY = '0.2'


class Genome:

    '''
    a reference genome from which base sequences can be read

    currently, the file that represents the reference genome can either
    be a simple sequence of the bases ``CAGT`` (``.bases`` file format)
    or a FASTA file. using a ``.bases`` file is slightly more efficient
    since the bases are not read into memory.
    '''

    def __init__(self, path, identifier=None, description=None):
        '''
        :param path: name of file to read bases from; can be ``.bases``
            file that directly contains base sequence (without any
            whitespace) or a file in FASTA format (only first genome is read)
        :param identifier: short identifier of genome; will be read from
            FASTA file if none specified
        :param description: text description; will also be read from
            FASTA file if none specified
        '''
        self.path = path
        self.f = file(path,'r')

        if self.f.read(1) == '>':
            self.fasta = True
            self.f.seek(0)
            defline = self.f.readline()
            idx = defline.find(' ')
            if identifier is None:
                if idx == -1:
                    identifier = defline[1:]
                else:
                    identifier = defline[1:idx]
            if description is None:
                if idx != -1 and idx < len(defline):
                    description = defline[idx + 1:]

            # read whole sequence into memory
            self.bases = ''
            self.bases = ''.join([line.rstrip('\n\r')
                    for line in self.f.readlines()])
            if '>' in self.bases:
                lo.info('%s contains several genomes; only first read' % path)
                self.bases = self.bases[:self.bases.index('>')]
            self.size = len(self.bases)
            self.f.close()
            lo.debug('read %d bytes FASTA sequence "%s" into memory' % (
                self.size, identifier))

        else:
            self.fasta = False
            self.f.seek(0,2)
            self.size = self.f.tell()

        self.identifier = identifier
        self.description = description

    def read(self, pos, length):
        '''
        :param pos: index of first base to read (starting at ``1``!)
        :param length: number of bases to read
        :returns: a base string of length ``length``
        '''
        if self.fasta:
            return self.bases[pos-1:pos-1 + length]
        self.f.seek(pos-1) # pos starts at 1...
        return self.f.read(length)

    def seq(self, start, stop, left=0, right=0, **kwargs):
        '''
        :param pos: index of first base to read (starting at ``1``!)
        :param stop: index of last base to read
        :param left: length of flank on the left
        :param right: length of flank on the right
        :param kwargs: additional parameters passed to :py:class:`.Sequence`
        :returns: a :py:class:`.Sequence` with length
            ``stop - start + 1 + left + right`` and ``left=left`` and
            ``right=right``
        '''
        bases = self.read(start-left, stop-start+1+left+right)
        return Sequence(bases, left, right, pos=start-left, **kwargs)

    def __str__(self):
        return self.identifier


class Gene:

    ''' defines a gene within a :py:class:`.Genome` '''

    def __init__(self, genome, identifier,
            start, stop, promoter_end=None, plus_strand=True,
            coding=True):
        '''
        :param genome: a :py:class:`.Genome`
        :param identifier: a short identifier of the gene
        :param start: position of the first coding base
        :param stop: position of the last coding base
        :param promoter_end: position of the base after the last
            base of the the promoter (defaults to ``start``)
        :param plus_strand: whether the gene is read from the ``+`` strand
        :param coding: whether the gene is coding (if set to ``False``, the
            method :py:meth:`mut2str` will report base changes instead of
            amino-acid changes
        '''
        self.genome = genome
        self.identifier = identifier
        self.plus_strand = plus_strand
        self.coding = coding
        assert start<=stop, 'start position must be smaller than stop position'
        #assert (stop-start +1)%3 == 0, 'gene length must be multiple of 3'
        self.start = start
        self.stop = stop
        if promoter_end is None:
            promoter_end = start
        self.promoter_end = promoter_end

    #TODO implement 2/3 mutations on same codon
    def mut2str(self, pos, newbase):
        ''' :param pos: refers to the absolute base position (within the genome)
            :param newbase: refers to the **plus strand**
            
            :returns: a string in the form ``'gene.XposY'`` where pos
                denominates the amino acid position and X is the symbol for
                wildtype and Y the symbol for the mutant amino acid;
                if the mutation is *before* the gene, ``'gene promoter mutation -X'``
                is returned; if the mutation is *after* the gene, ``'?'`` is
                returned; if the gene is *not coding*, then a string in the
                form 'posXY' is returned, where pos is the base position
                relative to the beginning of the gene, X will be the original base,
                and Y the new base '''

        if pos < self.promoter_end:
            return '%s promoter mutation %d' % (self.identifier, pos - self.promoter_end)
        elif pos < self.start or pos > self.stop:
            return '?'
            #return '??? (%d not within %d-%d)'%( pos, self.start, self.stop)

        pos1 = pos - self.start + 1
        codon_nr = (pos-self.start)/3 + 1
        codon_start = self.start+(codon_nr-1)*3
        codon_mut = pos-codon_start
        codon = self.genome.seq(codon_start, codon_start+2)
        oldbase = self.genome.read(pos, 1)

        if not self.plus_strand:
            pos1 = self.stop - pos + 1
            codon_nr = (self.stop-pos)/3 +1
            codon_mut = 2-codon_mut
            codon = codon.reverse()
            newbase = codon.pairs[newbase]
            oldbase = codon.pairs[oldbase]

        if self.coding:
            aa1 = codon.transcribe()
            aa2 = codon.transcribe(mutations=((codon_mut, newbase),))
            return self.identifier + '.' + aa1 + str(codon_nr) + aa2
        else:
            return self.identifier + '.' + str(pos1) + oldbase + newbase

    def __str__(self):
        ''' :returns: a string representation of the gene '''
        if self.plus_strand:
            return 'gene %s %d..%d'%(
                    self.identifier, self.start, self.stop)
        else:
            return 'gene %s complement(%d..%d)'%(
                    self.identifier, self.start, self.stop)


class Sequence(object):

    ''' a sequence of bases with a margin (``left`` and ``right``); the
        optional attribute ``pos`` refers to the **first base after
        the left flank** (i.e. the first base of the region of interest)

        when accessing the bases of the sequence, indexing starts with the
        **first base of the margin** (this corresponds to base position 
        ``pos-left`` within the genome) '''

    ''' how to map matching bases -- ``N`` can appear in some ``.fastq`` files '''
    pairs = { 'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N' }

    ''' DNA amino acid codon table, copied from
        https://en.wikipedia.org/wiki/DNA_codon_table '''
    code = { 
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
            'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'TGT':'C', 'TGC':'C', 'TGA':'$', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
            'TAT':'Y', 'TAC':'Y', 'TAA':'$', 'TAG':'$', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
            'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        }


    def __init__(self, bases, left=0, right=0, pos=None, plus_strand=True):
        ''' 
        :param bases: string of base letters ``ACGT``
        :param left: how many of theses bases are spacing on the left side
            (at 5' end = towards lower pos on plus strand)
        :param right: how many of theses bases are spacing on the right side
            (at 3' end = towards lower pos on minus strand)
        :param pos: position of the actual sequence (without flanks) within
            the gnome; note that ``bases[0]`` is positioned at ``pos-left``
        :param plus_strand: whether this sequence is on the plus strand

        initialize a :py:class:`Sequence` with a base sequence
        '''

        self.bases = bases
        self.left = left
        self.right = right
        self.pos = pos
        self.plus_strand = plus_strand

    def __len__(self):
        ''' :returns: the length of the sequence, //including// margins '''
        return len(self.bases)

    def __getitem__(self, idx):
        '''
        :returns: base at position ``idx``, starting from first base of
            left flank '''
        return self.bases[idx]

    def __setitem__(self, idx, value):
        ''' change single base in sequence; ``pos=0`` refers to first base
            of left flank '''
        if idx<0 or idx>len(self.bases)-1:
            raise IndexError
        self.bases = self.bases[:idx] + value + self.bases[idx+1:]

    def reverse(self):
        ''' :returns: the complementary sequence '''
        return Sequence(''.join([self.pairs[b]
                        for b in self.bases])[::-1],
                pos=self.pos, plus_strand=not self.plus_strand,
                left=self.left, right=self.right)

    def plus_idx(self, idx):
        ''' :returns: the index that corresponds to ``idx`` on the ``+`` strand '''
        if self.plus_strand:
            return idx
        else:
            return len(self.bases)-idx-1

    def plus_base(self, base):
        ''' :returns: the base that corresponds to ``base`` on the ``+`` strand '''
        if self.plus_strand:
            return base
        else:
            return self.pairs[base]

    def get_aa(self, pos, mutations=list()):
        '''
        :param pos: base position (with ``0`` corresponding to first base of left flank)
        :param mutations: list of ``[pos, base]`` of mutations to apply
        :returns: resulting amino acid at specified position
        '''
        pos0 = pos - pos%3
        codon = [c for c in self[pos0:pos0+3]]
        for pos, newbase in mutations:
            if pos-pos0<3 and pos-pos0>=0:
                codon[pos-pos0] = newbase
        return self.code[''.join(codon)]

    def transcribe(self, mutations=list()):
        '''
        :param mutations: list of ``[pos, base]`` of mutations to apply
        :returns: string of one-letter amino acid abbreviations corresponding to
            transcribed sequence after applying mutations
        '''
        ret = []
        for pos in range(len(self)/3):
            ret.append(self.get_aa(pos*3, mutations))
        return ''.join(ret)

    def apply_mutations(self, mutations):
        '''
        :param mutations: list of ``[pos, base]`` of mutations to apply

        applies mutations to ``.bases``
        '''
        bases = list(self.bases)
        for pos, newbase in mutations:
            bases[pos] = newbase
        self.bases = ''.join(bases)


class Template(object):

    ''' object with identifier that produces :py:class:`Sequence` '''

    def __init__(self, identifier):
        ''' subclasses must assure that the identifier is **unique**, i.e.
            that the same identifier always represents the same template
            (and thereby base sequence) because downstream code relies on
            reconstructing sequences from templates that are identified by
            this string '''
        self.identifier = identifier

    def validate(self, coverage):
        '''
        :param coverage: :py:class:`kvarq.analyse.Coverage`
        :returns: whether template was "found" in ``.fastq`` file
        '''
        #TODO make this dependent on coverage etc
        return coverage.mean(include_margins=False) >= 2

    def seq(self):
        ''' generate a :py:class:`.Sequence` from template '''
        raise NotImplementedError

    def __str__(self):
        ''' :returns: unique identifier '''
        return self.identifier


class StaticTemplate(Template):

    ''' a subclass of :py:class:`.Template` that cannot generate flanks '''

    def __init__(self, bases, identifier=None):
        ''' :param bases: base sequence of template '''
        if not identifier:
            identifier = bases
        super(StaticTemplate, self).__init__(identifier)
        self.bases = bases
        self.identifier = identifier

    def seq(self):
        return Sequence(self.bases)


class DynamicTemplate(Template):

    ''' a subclass of :py:class:`.Template` that can generate flanks '''

    def seq(self, spacing=0):
        '''
        :param spacing: length of flanks
        :returns: :py:class:`.Sequence` as read from the ``+`` strand
        '''
        raise NotImplementedError


class TemplateFromGenome(DynamicTemplate):

    '''
    generate a template from a reference genome

    because the reference genome is normally the sequence from a common ancestor,
    templates generated from it should be found in ``.fastq`` files unless a
    deletion occurred (contrary to :py:class:`.SNP`)
    '''

    def __init__(self, genome, start, stop, direction='+', aa_pos0=1, poslist=None):
        '''
        :param genome: a :py:class:`Genome` instance
        :param start: first base of sequence (start counting at 1)
        :param stop: last base of sequence (start counting at 1)
        :param direction: specifies from which strand the sequence is
            transcribed (either ``'+'`` or ``'-'``)
        :param aa_pos0: number of first amino acid transcribed (would be ``1``
            if template starts at beginning of gene)
        :param poslist: a sequence of positions that can later on be
            used e.g. for filtering of mutations (if only specific
            mutations are validated by the specified ``Test`` source)

        create a template from a genome
        '''

        identifier = '%s[%d:%d](%s)'%(str(genome), start, stop, direction)
        super(TemplateFromGenome, self).__init__(identifier)

        assert start<=stop
        assert direction in '+-'

        self.genome = genome
        self.start = start
        self.stop = stop
        self.aa_pos0 = aa_pos0
        self.direction = direction
        self.poslist = poslist

    def seq(self, spacing=0):
        return self.genome.seq(self.start, self.stop, spacing, spacing)

    #TODO? move into Gene
    def transcribe(self, mutations=None):
        ''' transcribes the sequence (from the strand specified by
            ``.direction`` '''
        seq = self.seq()
        if mutations:
            seq.apply_mutations(mutations)
        if self.direction == '-':
            seq = seq.reverse()
        return seq.transcribe()

    def mutations(self, coverage):
        ''' filters ``coverage.mutations`` and returns
            ``[[pos, base], ...]`` of the most prevalent mutation where
            ``pos`` is relative to ``.start`` '''
        ret = []

        mean = coverage.mean()
        std = coverage.std()

        for cpos, bases in coverage.mutations.items():
            # ignore mutations outside template region
            if cpos<coverage.start or cpos-coverage.start>=len(self.seq()):
                continue

            # pick most prevalent mutation
            basecounts = [
                    (base, sum([1 for b in bases if b==base]))
                    for base in set(bases)
                ]
            base, n = sorted(basecounts, key=lambda x:-x[1])[0]

            #TODO make this dependent of .fastq coverage etc
            if n>1 and n>mean-1.5*std:
                ret.append((cpos - coverage.start, base))

        return ret

    #TODO move into Gene
    def aa_mutations(self, mutations):
        ''' returns a list of ``[[aa_pos, old_aa, new_aa], ...]]`` from
            specified mutations ``[[pos, base], ...`` '''
        aa1 = self.transcribe()
        aa2 = self.transcribe(mutations)

        ret = []
        for i, old_aa in enumerate(aa1):
            if aa2[i] != old_aa:
                ret.append(( i+self.aa_pos0, old_aa, aa2[i] ))
        return ret


class SNP(TemplateFromGenome):

    '''
    represents a single nucleotide polymorphism

    note that, contrary to other :py:class:`.TemplateFromGenome`, the
    sequence of a SNP is the **mutant** version, i.e. if the template
    is not found in the ``.fastq`` file then the mutant defined by the
    SNP is not present (it might either be the wild type or another
    mutant)
    '''

    def __init__(self, genome, pos, base, orig=None, force=False):
        '''
        :param genome: reference genome from which to read bases
        :param pos: position of SNP within reference genome
        :param base: mutant base
        :param orig: original base
        :param force: do not check whether base/orig are compatible with
            bases read from genome
        '''

        super(SNP, self).__init__(genome, pos, pos)

        self.base = base
        self.orig = orig
        oldbase = self.genome.read(pos, 1)
        if not force:
            if orig: assert oldbase == self.orig, \
                    'expected orig %s found %s' % (self.orig, oldbase)
            assert base != oldbase
        self.identifier = 'SNP%d%s%s'%(pos,oldbase,base)

    def seq(self, spacing=0):
        seq = super(SNP, self).seq(spacing=spacing)
        seq[spacing] = self.base
        return seq

    def validate(self, coverage):
        ''' :returns: ``True`` if SNP is present, given ``coverage`` '''
        c = coverage.coverage[coverage.start]
        m = len(coverage.mutations.get(coverage.start, []))
        # TODO make this dependent of .fastq coverage etc
        return c>=2 and m<c/2


class Reference:

    '''
    represents a (litterature) reference where some kind of genetic information is defined
    '''

    def __init__(self, descr):
        self.descr = descr


class Genotype(object):

    def __init__(self, identifier, gene=None):
        '''
        :param identifier: identifier of this genotype
        :param gene: :py:class:`.Gene` linked to this genotype
        '''
        self.identifier = identifier
        self.gene = gene

    def __str__(self):
        return str(self.identifier)

    def __repr__(self):
        return '<%s : "%s">'%(self.__class__.__name__, self.identifier)


class Test(object):

    '''
    a test object links a :py:class:`.Template` to a :py:class:`.Genotype`

    the template will be used to scan the ``.fastq`` file and the genotype
    will be used to synthesize a result. see also :py:class:`.Testsuite`
    '''

    def __init__(self, template, genotype, reference):
        assert '::' not in str(template)
        self.template = template
        self.genotype = genotype
        self.reference = reference

    def __str__(self):
        return '%s::%s'%(self.genotype, self.template)


class AnalysisException(RuntimeError):
    ''' risen if error occurs during :py:meth:`kvarq.genes.Testsuite.analyse` '''

class Testsuite(object):

    '''
    interpretes features from a ``.fastq`` file using an array of :py:class:`.Test`

    this is a generic implementation that simply shows whether specified SNPs
    have occurred and whether there were any mutations found in the regions.
    subclassing testsuites its :py:meth:`._analyse` method in an application
    specific manner.
    '''

    def __init__(self, tests, version):
        '''
        :param tests: list of :py:class:`.Test` that is used by this testclass
        :param version: string used to make sure version loaded from ``.json``
            is compatible -- the minor number should be changed every time a
            change is made to the testsuite and the major number should be
            changed every time a non-backward-compatible change is made
            (such as the addition of a new test)
        '''
        self.tests = tests
        self.version = version

    def _analyse(self, coverages):
        ''' method doing the actual work for :py:meth:`analyse` to be overwritten
            by more intelligent subclasses '''
        ret = []

        for test in self.tests:
            coverage = coverages[test]
            seq = test.template.seq()

            if isinstance(test.template, SNP):
                if test.template.validate(coverage):
                    ret.append(str(test))

            elif isinstance(test.template, TemplateFromGenome):
                for pos, newbase in test.template.mutations(coverage):
                    oldbase = seq[pos]
                    ret.append('%d%s%s'%(
                            pos + test.template.start, oldbase, newbase))
                    if test.genotype.gene:
                        ret[-1] += '=' + test.genotype.gene.mut2str(
                                pos + test.template.start, newbase)
        return ret

    def analyse(self, analyser):
        '''
        :param analyser: :py:class:`kvarq.analyse.Analyser` containing the
            information gathered during a scanning process
        :returns: findings (can be string or array of strings)

        This method raises an :py:class:`.AnalysisException` if a
        :py:class:`.Test` that is included in this testsuite is not found in
        the ``analyser`` -- this could happen when a new ``Test`` is added to
        the testsuite without incrementing the ``TESTSUITE_VERSION`` minor number
        and then a previously saved ``.json`` file is decoded using the new
        version of the testsuite.

        No exception is risen if the analyser defines additional tests that
        were defined by a previous testsuite but are no longer needed by the
        current testsuite.

        .. automethod:: _analyse
        '''
        try:
            coverages = dict([(test, analyser[test]) for test in self.tests])
        except KeyError, e:
            raise AnalysisException('template "%s" not found' % str(test.template))
        return self._analyse(coverages)

    def __str__(self):
        return 'generic Testsuite with %d tests' % len(self.tests)


class TestsuiteLoadingException(Exception):
    ''' exception risen if error is encountered while loading a testsuite '''

def load_testsuite(fname):
    '''
    :param fname: path of ``.py`` testsuite file

    loads a modular testsuite from a file; see :ref:`testsuites`

    **beware** that the testsuite is a python file and can execute arbitrary
    code

    raises :py:class:`TestsuiteLoadingException` if file format of specified
    testsuite is invalid
    '''

    name = os.path.splitext(os.path.basename(fname))[0]
    if '-' in name:
        name = name[:name.index('-')]
    namespace = dict(
            __file__=fname,
            __module__='kvarq.testsuites.' + name
        )

    try:
        sys.path.insert(0, os.path.dirname(fname))
        execfile(fname, namespace)
        del sys.path[0]
    except Exception as e:
        raise TestsuiteLoadingException('exception while reading file : %s [%s]' % (
            str(e), format_traceback(sys.exc_info())))

    if not 'GENES_COMPATIBILITY' in namespace:
        raise TestsuiteLoadingException('module defines no "GENES_COMPATIBILITY"')
    
    compat = StrictVersion(namespace['GENES_COMPATIBILITY'])
    version = StrictVersion(COMPATIBILITY)

    if compat > version or compat.version[0] != version.version[0]:
        raise TestsuiteLoadingException('incompatible : %s needed, got %s' %
                (compat, version))

    if not name in namespace:
        raise TestsuiteLoadingException('could not import "%s"' % name)
    if not isinstance(namespace[name], Testsuite):
        raise TestsuiteLoadingException('modules defines "%s" but is of type %s' %
                type(namespace[name]))

    return namespace[name]

