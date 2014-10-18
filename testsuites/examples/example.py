
# this is an example testsuite that illustrates how to write simple
# SNP/region based testsuite to be used with kvarq

# the testsuite can be included during the scanning by using the
# command line parameter '-t' or in the configuration window in the GUI

# see the kvarq documentation for more information:
# http://kvarq.readthedocs.org/en/latest/testsuites.html


# the version specifies the version of the testsuite itself; this version
# string is included in the .json scan results
# the minor number should be increased every time the file changes. the major
# number should be increased when the changes are not backwards compatible
# (e.g. when a new test is added)

VERSION = '0.1'

# this version is compared against the COMPATIBILITY global module variable
# defined in kvarq.genes
# as before, compatibility is warranted if the first number is equal and the
# second equal or lower (than the one defined in kvarq.genes)

GENES_COMPATIBILITY = '0.0'


# we use these classes to define our testsuite
from kvarq.genes import Genotype, Gene, Test, Testsuite, Reference, SNP, TemplateFromGenome

# load hypothetical MTB ancestor genome from '../MTBC' directory
# (shipped together with KvarQ)
from kvarq.genes import Genome
import os.path
MTBC_dir = os.path.join(os.path.dirname(__file__), os.pardir, 'MTBC')
ancestor_path = os.path.join(MTBC_dir, 'MTB_ancestor_reference.bases')
ancestor = Genome(ancestor_path, 'MTB ancestor')

# use this for loggging (displayed on console / in main GUI window)
from kvarq.log import lo


# references tell where more information ont he mutations can be found
tbdream = Reference('TBDReamDB : see http://tbdreamdb.com/')

# the first genotype simply signals isoniazid resistance
inhA = Genotype('Isoniazid resistance')
# the second genotype also signals isoniazid resistance but indicates
# the gene to which it belongs to -- this enables output of resistance
# mutation in the familiar gene.XposX format
katG = Genotype('Isoniazid resistance', Gene(ancestor,'katG', 2153889, 2156111, plus_strand=False))

# define two SNPs : 1673432TA and 1673432TC -- only specified mutations will
# be found (i.e. 1673432TATG would not be reported)
# note that the SNP is simply the "template" for that will be used when scanning
# for mutations in the FastQ reads; the "test" as a whole defines a template,
# a genotype (inhA in this case) and the resource from the information is drawn
SNP1 = Test(SNP(genome=ancestor, pos=1673432, orig='T', base='A'), inhA, tbdream)
SNP2 = Test(SNP(genome=ancestor, pos=1673432, orig='T', base='C'), inhA, tbdream)

# define a (short) region that should be scanned for ANY mutations here we're
# interested in the codon 2155167-2155169; by specifying where the gene is read
# from (minus strand) and the position of the amino acid is produced by this
# codon (in this case the gene starts at 2153889, therefore the amino acid is
# ((2155167-2153889)/3 +1)=427) it is later possible to check for (non)
# synonymous mutations as before, the "test" consists of a template, a genotype
# and a resource (but the "template" is a region and not a SNP as before)
katG_codon = Test(TemplateFromGenome(genome=ancestor, start=2155167,
    stop=2155169, direction='-', aa_pos0=(2155167-2153889)/3 +1), katG,
    tbdream)


# it's important to NAME the testsuite the SAME AS THE FILENAME up to the first
# dash !  (e.g. it's possible to rename this file to "example-0.1.py")
example = Testsuite([SNP1, SNP2, katG_codon], VERSION)


# note that this testsuite is very simple and will simply eport any mutations
# found in the FastQ file -- often it makes sense to subclass the Testsuite
# class (and redefine the _analyse method) to get a fine-grained control on how
# the mutations are synthesized into a result... see the source code of
# kvarq.genes.phylo as an example

