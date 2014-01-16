
from kvarq import genes
from kvarq.genes import TemplateFromGenome, Gene, load_testsuite, Genome
from kvarq import analyse

import unittest
import os.path

MTBCpath = os.path.join(os.path.dirname(__file__), os.path.pardir, 'testsuites', 'MTBC')
ancestor = Genome(os.path.join(MTBCpath, 'MTB_ancestor_reference.bases'), 'MTB ancestor')

class GenesTest(unittest.TestCase):


    def test_gene(self):
        rpoB = Gene(genome=ancestor, identifier='rpoB',
                start=759807, stop=763325, plus_strand=True)
        assert rpoB.mut2str(761082, 'C') == 'rpoB.G426R'

        MTB10 = Gene(genome=ancestor, identifier='MTB10',
                start=1021344, stop=1021643, plus_strand=False)
        assert MTB10.mut2str(1021600, 'C') == 'MTB10.T15R'
        assert MTB10.mut2str(1021601, 'A') == 'MTB10.T15S'


    def test_SNP(self):

        x = genes.SNP(genome=ancestor, pos=3920109, base='T')
        assert x.seq(spacing=3).bases  == 'CGATATT'

        x = genes.SNP(genome=ancestor, pos=3920109, base='T')
        assert x.seq(spacing=12).bases == 'TTACTGCGCCGATATTCGCACACCT'

        x = genes.SNP(genome=ancestor, pos=2427828, base='C')
        assert x.seq(spacing=12).bases == 'CCACAGTGTGAGCCCTAGTCCGACG'


    def test_reverse(self):
        assert genes.Sequence('AAACGT').reverse().bases == 'ACGTTT'


    def test_code(self):
        seq = genes.Sequence('GCTTGTGATTGC')
        for i in range(4):
            for j in range(3):
                assert seq.get_aa(i*3+j) == 'ACDC'[i]
        assert seq.get_aa(1, [(1, 'T')]) == 'V'
        assert seq.transcribe() == 'ACDC'

        # test forward coding template
        Rv0880 = TemplateFromGenome(ancestor, 978934, 979365, direction='+')
        Rv0880seq = Rv0880.seq()
        assert Rv0880seq.plus_strand
        assert Rv0880seq.bases.startswith('GTGCTTGACAGCGA')
        assert Rv0880seq.transcribe().startswith('VLDSDARLASDL')
        assert Rv0880seq.transcribe(mutations=((1, 'G'),)).startswith('GLDSDARLASDL')

        # test reverse coding template
        Rv0883c = TemplateFromGenome(ancestor, 980506, 981267, direction='-')
        Rv0883cseq = Rv0883c.seq() # returns sequence from '+' strand !
        assert Rv0883cseq.plus_strand
        assert Rv0883cseq.bases.startswith('CTAGCGACG')
        assert Rv0883c.transcribe().startswith('MRELKVVGLD')
        pos = len(Rv0883cseq) -2 -1 # second last base
        assert Rv0883c.transcribe(mutations=((pos, 'G'),)).startswith('IRELKVVGLD')


    def test_mutations(self):
        # first test mutation validation of SNP
        snp1000 = genes.SNP(ancestor, 1000, base='C', orig='G')
        seq = snp1000.seq(spacing=25)
        coverage = analyse.Coverage(seq)
        # no mutations
        assert not snp1000.validate(coverage)
        # not enough coverage
        coverage.mutations = dict([(25, 'C')])
        assert not snp1000.validate(coverage)
        # not enough mutations
        coverage.coverage = [20] * len(coverage.coverage)
        coverage.mutations = dict([(25, 'C'*10)])
        assert not snp1000.validate(coverage)
        # this should validate
        coverage.mutations = dict([(25, 'C')])
        assert snp1000.validate(coverage)

        # then test mutation filtering of TemplateFromGenome
        embB = genes.TemplateFromGenome(ancestor, 4246514, 4249810, direction='+')
        seq = embB.seq(spacing=25)
        coverage = analyse.Coverage(seq)
        coverage.mutations = dict([
                (25 + 0, 'TGC'), # 'A' from 'ATG'=M (codon=1) : no mutation
                (25 + 4, 'GGGGGGGGGGGAT'), # 'A' from 'ACA'=T (codon=2) : -> 'AGA'=R
            ])
        # mutation not validated when coverage very high
        coverage.coverage = [1000] * len(coverage.coverage)
        assert len(embB.mutations(coverage)) == 0
        # one mutation detected with coverage == mutation count
        coverage.coverage = [10] * len(coverage.coverage)
        mutations = embB.mutations(coverage)
        assert len(mutations) == 1
        assert mutations[0] == (4, 'G')
        # this should be translated into "T2R"
        aa_mutations = embB.aa_mutations(mutations)
        assert len(aa_mutations) == 1
        assert aa_mutations[0] == (2, 'T', 'R')

if __name__ == '__main__': unittest.main()

