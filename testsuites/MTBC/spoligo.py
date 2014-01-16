'''
defines spoligo sequences (taken from SpolPred by Francesc Coll) and some helper
functions to handle them
'''

VERSION = '0.2'
from kvarq.genes import COMPATIBILITY as GENES_COMPATIBILITY

from kvarq.genes import Genotype, Test, StaticTemplate, Reference, Testsuite


def code(spnrs):
    ''' spoligo0-spoligo42 can be represented as a 15 digit code
    
        the first 14 digits are octal representation of 3 spoligos each (a 4 as
        the first digit is equal to spoligo0+, spoligo1-, spoligo2-) -- the last
        digit is binary for spoligo42 '''

    value = 0
    spoligo42 = '0'
    for spnr in spnrs:
        if spnr == 42:
            spoligo42 = '1'
        else:
            value += 2**(41-spnr)
    octal = oct(value)[1:].rstrip('L')
    octal = '0'*(14-len(octal)) + octal

    return octal + spoligo42


class SpoligoTestsuite(Testsuite):

    def __str__(self):
        return 'TB spoligos'

    def _analyse(self, coverages):

        spnrs = [spnr for spnr, spoligo in enumerate(self.tests)
                if spoligo.template.validate(coverages[spoligo])]

        if not spnrs or sum([coverages[self.tests[spnr]].mean() for snpr in spnrs]) / len(spnrs) < 10:
            remark = ' -- low coverage (mean below 10x)'
        else:
            remark = ''

        spoct = code(spnrs)

        spbin = ''.join(
                [bin(int(x))[2:].rjust(3, '0') for x in spoct[:14]] +
                [bin(int(x))[2:] for x in spoct[14:]]
            )
        spnum = ','.join([str(spnr) for spnr in sorted(spnrs)])

        return ' '.join([spoct, spbin]) + remark


class Spoligo(Genotype):

    def __init__(self, number):
        assert number >= 0 and number <= 42
        super(Spoligo, self).__init__('spoligo'+str(number))
        self.number = number


spolpred = Reference('SpolPred.cpp by Francesc Coll')

spoligo_counter = 0
def NextSpoligo():
    global spoligo_counter
    ret = Spoligo(spoligo_counter)
    spoligo_counter += 1
    return ret

spoligo_seqs = [
    Test(StaticTemplate('TGATCCAGAGCCGGCGACCCTCTAT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CAAAAGCTGTCGCCCAAGCATGAGG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TAGAAGGCGATCACTGGAAGCACGG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CTGATGATTGGTCGGCGTATGACGT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TAATCCCGCACAAGTGGTCAGAAAA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GAAATTGAAGCCGGAAATGACGACG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GCAGCCCCGAGTACTCGCTCTCCTC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CGGCGAGGCTGGGGGCGGTTTCACG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GCTGTCAGCACATGGGATTCCGAGT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GGAAGTCAACTAGAGCGGGTGTCGA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CCAGGTTGCCGCCGCCGTTGCTCAC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('ATCTCCCCGGGCGGGCAGCAGATAT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GGGAGAGGGAATGGCAATGATGGTC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CCGAGCCGACCATCCGCATCACACC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CGAAATTCACTGCGCGTTATTCAAG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GATTTACGACGCTGACGGGAACTCG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CGGAGTCATCCGCGCGGGCCGGCGC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CATCTGCAGCTCGCCCGGGTCCATG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('ACCAGGATCAGCGCCAAGCCAGTTA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TGATCTTCTCTCCTGGCGAGGTCAA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TCGACGATTGGGACATCGACATCGA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TTGTCTCAATCGTGCCGTCTGCGGT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CGAGCTGGACCGCATCAGCGATGCT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CGAGCACGTCTCACCCAGCAGGCGG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TGACAGGGTGCGGTGGTCGCTGATC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GCGCCGGATGATGGTGGTGCTGAAG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('ATCCGCGGGAAGAGATCACGAATCC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GTTGTGATCGCTAAACGCCGGGGCA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TGGTCGTGTCGTGGAGCCTGTATTT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GGCTGGAAAAGGGCGCGGGGCAACC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('ACTTGATCGACGCGAACCTGTCTGA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TGAACACGCCGATACCTATTTGGTC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TCAAGTGCGGCACCGCCGTCATGTC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TTCGACGGTGTGGGCGAGGTGACTT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GTTGGAAGCGTTTCGAGCGTACGGA'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GCTGCGGATGTGGTGCTGGATTTCG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('AAGGGGGACTGTGGACGAGTTCGCG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GCGCACAACGCATCCGCCATCCACG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CCACGCCGATTTACTGGCCATCGTC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GGACCTGTATGAGGCACAGATGGCG'), NextSpoligo(), spolpred),
    Test(StaticTemplate('TACCTGATAGAAGCCGGAAAGCTCC'), NextSpoligo(), spolpred),
    Test(StaticTemplate('GTCGCGCTCGTCCATGTCCCACCAT'), NextSpoligo(), spolpred),
    Test(StaticTemplate('CTCCCGCACCCGGTGCGATTCTGCG'), NextSpoligo(), spolpred),
]

assert len(spoligo_seqs) == 43

spoligo = SpoligoTestsuite(spoligo_seqs, VERSION)
