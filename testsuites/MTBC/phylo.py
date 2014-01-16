
VERSION = '0.6'
from kvarq.genes import COMPATIBILITY as GENES_COMPATIBILITY

from kvarq.log import lo
from kvarq.genes import Reference, Test, SNP, Genotype, Testsuite

from _util import ancestor

class PhyloTestsuite(Testsuite):

    def __str__(self):
        return 'TB lineage SNPs'

    def score_SNPs(self, genotypes, coverages):
        ''' returns a dictionary ``{genotype:scores, ...}`` where ``scores`` is an
            array where every SNP for the given genotype is represented with a
            boolean value indicating whether the SNP is covered "sufficiently" 
            to be considered positive '''

        ret = {}

        for test in self.tests:

            coverage = coverages[test]
            genotype = test.genotype

            if genotype in genotypes:
                if test.template.validate(coverage):
                    ret.setdefault(genotype, []).append(True)
                else:
                    ret.setdefault(genotype, []).append(False)

        return ret

    def _analyse(self, coverages):
        mls = []

        #TODO choose criteria dynamically

        for ml, xs in self.score_SNPs(Lineage.roots, coverages).items():
            lo.debug(str(ml)+' : '+str(xs))

            if sum(xs)>1:
                # we need at least two positive SNPs
                mls.append(ml.name)

                if 0 in xs:
                    # flag if one of the SNPs is not found
                    # mls[-1] += ' (?)'
                    pass

                if ml.children:
                    sls = []

                    # co-mutants complicate our life somewhat
                    slsc = self.score_SNPs(ml.children, coverages)
                    slsc_byname = {}
                    slsc_comutants = {}
                    for sl, xs_ in slsc.items():
                        slsc_byname.setdefault(sl.name, []).extend(xs_)
                        if sl.comutant:
                            slsc_comutants.setdefault(sl.name, []).extend([sl.comutant] * sum(xs_))

                    for slname, xs_ in slsc_byname.items():
                        comutants = ''.join(slsc_comutants.get(slname, []))
                        lo.debug('sublineage '+slname+' : '+str(xs_)+' comutants '+comutants)
                        if sum(xs_)>1:
                            sls.append(slname)
                            if comutants:
                                sls[-1] += '_' + comutants

#                            if 0 in xs_: # does not make sense when using comutants
#                                sls[-1] += ' (?)'

                    if sls:
                        mls[-1] += '/' + '-'.join(sls)

        coverages = sorted([coverage.mean(include_margins=False)
                for coverage in coverages.values()])
        remark = ''
        if coverages[len(coverages)/2] < 10:
            remark = ' -- low coverage (median below 10x)'

        if not mls:
            return '?' + remark

        return ' // '.join(mls) + remark


class Lineage(Genotype):

    roots = []

    def __init__(self, name, parent=None, color=None, origin=None, comutant=None):
        super(Lineage, self).__init__(name)
        self.name = name
        self.parent = parent
        self.color = color
        self.origin = origin
        self.comutant = comutant
        self.children = []

        if parent:
            parent.children.append(self)
            if color==None:
                color = self.parent.color
        else:
            Lineage.roots.append(self)


comas09  = Reference('PLoS ONE 2009 - Comas (monomorphic)')
stucki12 = Reference('Stucki et al. PLoS ONE 2012')

lineage1         = Lineage('lineage 1',color='magenta', origin='east africa, indian ocean, phillipines')
lineage2         = Lineage('lineage 2',color='blue', origin='TODO')
lineage_beijing  = Lineage('beijing sublineage',lineage2)
lineage3         = Lineage('lineage 3',color='purple', origin='east africa, central asia')
lineage4         = Lineage('lineage 4',color='red', origin='europe, america, africa')
lineage5         = Lineage('lineage 5',color='brown', origin='west africa 1')
lineage6         = Lineage('lineage 6',color='green', origin='west africa 2')
lineage7         = Lineage('lineage 7',color='yellow', origin='aethiopian')
lineage_animal   = Lineage('animal lineage')


mainlineages_SNPs = [
    Test(SNP(genome=ancestor,pos=3920109,base='T'),lineage1,stucki12),
    Test(SNP(genome=ancestor,pos=3597682,base='T'),lineage1,stucki12),
    Test(SNP(genome=ancestor,pos=1590555,base='T'),lineage1,stucki12),
    Test(SNP(genome=ancestor,pos=1834177,base='C'),lineage2,stucki12),
    Test(SNP(genome=ancestor,pos=3304966,base='A'),lineage2,stucki12),
    Test(SNP(genome=ancestor,pos=2711722,base='G'),lineage2,comas09),
    Test(SNP(genome=ancestor,pos= 301341,base='A'),lineage3,stucki12),
    Test(SNP(genome=ancestor,pos=4266647,base='G'),lineage3,stucki12),
    Test(SNP(genome=ancestor,pos= 157129,base='T'),lineage3,comas09),
    Test(SNP(genome=ancestor,pos=3326554,base='A'),lineage4,stucki12), # same as in ancestor
    Test(SNP(genome=ancestor,pos=2154724,base='C'),lineage4,stucki12), # same as in ancestor
    Test(SNP(genome=ancestor,pos= 648856,base='T'),lineage4,stucki12), # same as in ancestor
    Test(SNP(genome=ancestor,pos=1377185,base='G'),lineage5,stucki12),
    Test(SNP(genome=ancestor,pos= 801959,base='T'),lineage5,stucki12),
    Test(SNP(genome=ancestor,pos=2859147,base='T'),lineage5,stucki12),
    Test(SNP(genome=ancestor,pos=2427828,base='C'),lineage6,stucki12),
    Test(SNP(genome=ancestor,pos= 378404,base='A'),lineage6,stucki12),
    Test(SNP(genome=ancestor,pos=4269522,base='A'),lineage6,stucki12),
    Test(SNP(genome=ancestor,pos=  14806,base='C'),lineage7,stucki12),
    Test(SNP(genome=ancestor,pos=1663221,base='G'),lineage7,stucki12),
    Test(SNP(genome=ancestor,pos= 497126,base='A'),lineage7,stucki12),
    Test(SNP(genome=ancestor,pos=3480645,base='G'),lineage_animal,stucki12),
    Test(SNP(genome=ancestor,pos=1427476,base='T'),lineage_animal,stucki12),
    Test(SNP(genome=ancestor,pos=3624593,base='T'),lineage_animal,stucki12),
]

beijing_SNPs = [

    Test(SNP(genome=ancestor,pos=2112832,base='C'),lineage_beijing,stucki12),
    Test(SNP(genome=ancestor,pos=3587446,base='A'),lineage_beijing,stucki12),
    Test(SNP(genome=ancestor,pos=1849051,base='T'),lineage_beijing,stucki12),
]

phylo = PhyloTestsuite( mainlineages_SNPs + beijing_SNPs, VERSION )

