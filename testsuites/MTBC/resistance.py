
VERSION = '2.0'
from kvarq.genes import COMPATIBILITY as GENES_COMPATIBILITY

from kvarq.genes import Genotype, Test, Reference, SNP, Testsuite
from kvarq.genes import TemplateFromGenome, Gene

from _util import ancestor

class ResistanceTestsuite(Testsuite):

    def __str__(self):
        return 'TB resistance mutations & mutation regions'

    def _analyse(self, coverages):
        ret = []

        # will be set to True if any of the regions / SNPs has maximum
        # variant base below 90% (no matter whether templates validate or not)
        mixed = False

        for test in self.tests:

            coverage = coverages[test]
            seq = test.template.seq()

            # a) SNPs
            if isinstance(test.template, SNP):
                minf = coverage.minf()
                if minf < 0.9:
                    mixed = True
                if test.template.validate(coverage):
                    ret.append(str(test))
                    ret[-1] += '=' + test.genotype.gene.mut2str(
                            test.template.start, test.template.base)
                    # show percentage of most prominent base at SNP position
                    # if it's below 90%
                    if minf < 0.9:
                        ret[-1] += ' (' + str(int(100 * minf)) + '%)'
                continue

            # b) regions
            if not mixed and coverage.minf() < 0.9:
                mixed = True
            mutations = test.template.mutations(coverage)
            output = []
            for pos, newbase in mutations:

                oldbase = seq[pos]
                output.append('%d%s%s'%(
                        pos + test.template.start, oldbase, newbase))

                if test.genotype.gene:
                    output[-1] += '=' + test.genotype.gene.mut2str(
                            pos + test.template.start, newbase)

                mutf = coverage.fractions_at(pos).values()[0]
                if mutf < 0.9:
                    output[-1] += ' (' + str(int(100 * minf)) + '%)'

            aa1 = test.template.transcribe()
            aa2 = test.template.transcribe(mutations)

            # ignore resistance mutations without aa change
            if test.genotype.gene.coding and aa1 == aa2:
                continue

            # notify if mutation is not in "poslist"
            if hasattr(test.template, 'poslist'):
                poslist = test.template.poslist
                if poslist and not [m for m in mutations
                        if m[0]+test.template.start in poslist]:
                    output.append('[NONE OF MUTATIONS DOCUMENTED IN REFERENCE]')

            ret.append(str(test.genotype) + ' [' + ' '.join(output) + ']')

        RRDR_tests = [test for test in self.tests if test.genotype==RRDR]
        assert len(RRDR_tests) == 1
        if coverages[RRDR_tests[0]].mean(include_margins=False) < 10:
            ret.append('remark: low coverage (RRDR below 10x)')
        if mixed:
            ret.append('remark: mixed coverage')

        return ret


class DrugResistance(Genotype):

    def __init__(self, drug, gene, remarks=None):
        identifier = drug + ' resistance'
        if remarks:
            identifier += ' ('+remarks+')'
        super(DrugResistance, self).__init__(identifier)
        self.drug = drug
        self.gene = gene
        self.remarks = remarks


# MDR : rifampicin + any of (isoniazid, ?)
# XDR : MDR + fluoroquinolone + injectable (aminoglycosides)

inhA = DrugResistance('Isoniazid', Gene(ancestor,'inhA', 1674202, 1675011, promoter_end=1673440))
katG = DrugResistance('Isoniazid', Gene(ancestor,'katG', 2153889, 2156111, plus_strand=False))

RRDR = DrugResistance('Rifampicin', Gene(ancestor,'rpoB', 759807, 763325), 'RRDR')
rpoA = DrugResistance('Rifampicin', Gene(ancestor,'rpoA', 3877464, 3878507, plus_strand=False), 'compensatory')
rpoC = DrugResistance('Rifampicin', Gene(ancestor,'rpoC', 763370, 767320), 'compensatory')

QRDR = DrugResistance('Fluoroquinolones', Gene(ancestor,'gyrA', 7302, 9818), 'QRDR')
gyrA = DrugResistance('Fluoroquinolones', Gene(ancestor,'gyrA', 7302, 9818))
gyrB = DrugResistance('Fluoroquinolones', Gene(ancestor,'gyrB', 5123, 7267))

rpsL = DrugResistance('Streptomycin', Gene(ancestor,'rpsL', 781560, 781934))

#TODO report position only, is in ribosomal RNA ...
rrsS = DrugResistance('Streptomycin', Gene(ancestor,'rrsS', 1471846, 1473382))
rrsK = DrugResistance('Kanamycin/Amikacin', Gene(ancestor,'rrsK', 1471846, 1473382))

embB = DrugResistance('Ethambutol', Gene(ancestor,'embB', 4246514, 4249810))

pncA = DrugResistance('Pyrazinamide', Gene(ancestor,'pncA', 2288681, 2289241, plus_strand=False))


comas12 = Reference('Comas et al 2012 Nat Gen: Compensatory mutations...')
ramaswamy98 = Reference('Ramaswamy et al., Tuber Lung Dis 1998')
sun08 = Reference('Sun et al., Antimicr Agents 2008')
tbdream = Reference('TBDReamDB')
sebastien = Reference('Sebastien')
sebastien_= Reference('Sebastien ?')
david = Reference('David')



resistance_SNPs = [

    Test(SNP(genome=ancestor, pos=2155276, orig='C', base='T'), katG, tbdream),
    Test(SNP(genome=ancestor, pos=1673432, orig='T', base='A'), inhA, tbdream),
    Test(SNP(genome=ancestor, pos=1673432, orig='T', base='C'), inhA, tbdream),
    Test(SNP(genome=ancestor, pos=1673425, orig='C', base='T'), inhA, tbdream),
    Test(SNP(genome=ancestor, pos=3877949, orig='T', base='C'), rpoA, comas12),
    Test(SNP(genome=ancestor, pos=3877949, orig='T', base='G'), rpoA, comas12),
    Test(SNP(genome=ancestor, pos=3877960, orig='A', base='G'), rpoA, comas12),
    Test(SNP(genome=ancestor, pos=3877960, orig='A', base='C'), rpoA, comas12),
    Test(SNP(genome=ancestor, pos=764669, orig='C', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764670, orig='C', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764817, orig='T', base='C'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764817, orig='T', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764819, orig='T', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764822, orig='G', base='A'), rpoC, comas12), # changed orig T->G
    Test(SNP(genome=ancestor, pos=764822, orig='G', base='C'), rpoC, comas12), # changed orig T->G
    Test(SNP(genome=ancestor, pos=764840, orig='A', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764841, orig='T', base='C'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=764918, orig='G', base='C'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=765461, orig='A', base='C'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=765462, orig='A', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=765463, orig='C', base='G'), rpoC, comas12),
    Test(SNP(genome=ancestor, pos=7606, orig='C', base='A'), gyrA, tbdream),
    Test(SNP(genome=ancestor, pos=7677, orig='G', base='A'), gyrA, tbdream),
    Test(SNP(genome=ancestor, pos=7678, orig='C', base='G'), gyrA, tbdream),
    Test(SNP(genome=ancestor, pos=6767, orig='G', base='A'), gyrB, tbdream),
    Test(SNP(genome=ancestor, pos=6768, orig='G', base='A'), gyrB, tbdream),
    Test(SNP(genome=ancestor, pos=781687, orig='A', base='G'), rpsL, tbdream),
    Test(SNP(genome=ancestor, pos=781822, orig='A', base='C'), rpsL, tbdream),
    Test(SNP(genome=ancestor, pos=781822, orig='A', base='T'), rpsL, tbdream),
    Test(SNP(genome=ancestor, pos=781822, orig='A', base='G'), rpsL, tbdream),
    Test(SNP(genome=ancestor, pos=1472337, orig='C', base='A'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472337, orig='C', base='G'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472337, orig='C', base='T'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472358, orig='C', base='A'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472358, orig='C', base='G'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472358, orig='C', base='T'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472359, orig='A', base='C'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472359, orig='A', base='G'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472359, orig='A', base='T'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472362, orig='C', base='A'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472362, orig='C', base='G'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472362, orig='C', base='T'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472752, orig='A', base='C'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472752, orig='A', base='G'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1472752, orig='A', base='T'), rrsS, tbdream),
    Test(SNP(genome=ancestor, pos=1473246, orig='A', base='C'), rrsK, tbdream),
    Test(SNP(genome=ancestor, pos=1473246, orig='A', base='G'), rrsK, tbdream),
    Test(SNP(genome=ancestor, pos=1473246, orig='A', base='T'), rrsK, tbdream),
    Test(SNP(genome=ancestor, pos=1473247, orig='C', base='A'), rrsK, tbdream),
    Test(SNP(genome=ancestor, pos=1473247, orig='C', base='G'), rrsK, tbdream),
    Test(SNP(genome=ancestor, pos=1473247, orig='C', base='T'), rrsK, tbdream),
    Test(SNP(genome=ancestor, pos=4247429, orig='A', base='G'), embB, tbdream),
    Test(SNP(genome=ancestor, pos=4247431, orig='G', base='A'), embB, tbdream),
    Test(SNP(genome=ancestor, pos=4247431, orig='G', base='T'), embB, tbdream),
    Test(SNP(genome=ancestor, pos=4247431, orig='G', base='C'), embB, tbdream),
    Test(SNP(genome=ancestor, pos=4247429, orig='A', base='C'), embB, tbdream),
    Test(SNP(genome=ancestor, pos=4247730, orig='G', base='C'), embB, tbdream),
    Test(SNP(genome=ancestor, pos=4248003, orig='A', base='G'), embB, tbdream),
]

resistance_regions = [

    Test(TemplateFromGenome(genome=ancestor, start=2155167, stop=2155169, direction='-', aa_pos0=(2155167-2153889)/3 +1), katG, ramaswamy98),
    Test(TemplateFromGenome(genome=ancestor, start=761082, stop=761162), RRDR, ramaswamy98),
    Test(TemplateFromGenome(genome=ancestor, start=7521, stop=7583, poslist=[7521, 7522, 7523, 7569, 7570, 7571, 7572, 7573, 7574, 7581, 7582, 7583]), QRDR, sun08),
    Test(TemplateFromGenome(genome=ancestor, start=2288681, stop=2289241, direction='-'), pncA, david),

]


resistance = ResistanceTestsuite( resistance_SNPs + resistance_regions, VERSION )

