
VERSION = '0.2'
GENES_COMPATIBILITY = '0.1'

import os.path, functools

from kvarq.genes import Gene, Genotype, Test, Testsuite, Reference, SNP
from kvarq.genes import Genome

# old ebola genome from previous outbreak
EBOV76 = Genome(os.path.join(os.path.dirname(__file__), 'EBOV76.fasta'))


gire14 = Reference('Gire et al (2014) doi 10.1126/science.1259657')

# sub-lineages as defined in gire14
SL1 = Genotype('SL1')
SL2 = Genotype('SL2')
SL3 = Genotype('SL3')

# SNPs extracted from primary data using suppl/_extract_SNPs.py
SNPs = [
        Test(SNP(genome=EBOV76, pos=800, orig='C', base='T'), SL2, gire14),
        Test(SNP(genome=EBOV76, pos=1849, orig='T', base='C'), SL1, gire14),
        Test(SNP(genome=EBOV76, pos=6283, orig='C', base='T'), SL1, gire14),
        Test(SNP(genome=EBOV76, pos=8928, orig='A', base='C'), SL2, gire14),
        Test(SNP(genome=EBOV76, pos=10218, orig='G', base='A'), SL3, gire14),
        Test(SNP(genome=EBOV76, pos=13856, orig='A', base='G'), SL1, gire14),
        Test(SNP(genome=EBOV76, pos=15660, orig='T', base='C'), SL1, gire14),
        Test(SNP(genome=EBOV76, pos=15963, orig='G', base='A'), SL2, gire14),
        Test(SNP(genome=EBOV76, pos=17142, orig='T', base='C'), SL2, gire14),
    ]


class CountGenotype(Testsuite):

    def __str__(self):
        return 'Counting SNPs by genotype'

    def _analyse(self, coverages):

        found = dict()
        total = dict()

        for test in self.tests:
            if isinstance(test.template, SNP):
                ident = test.genotype.identifier
                total[ident] = total.get(ident, 0) + 1
                found[ident] = found.get(ident, 0) + (
                        test.template.validate(coverages[test]) and 1)

        return [
                '%s:%d/%d' % (ident, found[ident], total[ident])
                for ident in total
            ]

sierraleone14 = CountGenotype(SNPs, VERSION)

