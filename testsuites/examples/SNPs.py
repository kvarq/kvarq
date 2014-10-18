
# see ./example.py for explanation of general structure of testsuite

VERSION = '0.0'
GENES_COMPATIBILITY = '0.0'

import os.path

from kvarq.genes import Genome, Reference, SNP, Test, Testsuite, Genotype

def tsv2SNPs(path, genome, reference):

    tests = []
    for line in file(path):

        parts = line.strip().split('\t')
        name = parts[0]
        pos = int(parts[1])
        bases = parts[2].split('/')

        snp = SNP(genome=genome, pos=pos, orig=bases[0], base=bases[1])
        test = Test(snp, Genotype(name), reference)
        tests.append(test)

    return tests

here = os.path.dirname(__file__)
# we borrow the reference from the ../MTBC testsuites
genome_path = os.path.join(here, os.path.pardir, 'MTBC', 'MTB_ancestor_reference.bases')
genome = Genome(genome_path, 'MTB ancestor')
ref = Reference('specify reference here')
# load SNP information from .tsv file (can be edited with Excel)
SNPs = tsv2SNPs(os.path.join(here, 'SNPs.tsv'), genome, ref)

SNPs = Testsuite(SNPs, VERSION)

