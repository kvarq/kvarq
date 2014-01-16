
VERSION = '0.5'
GENES_COMPATIBILITY = '0.0'

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
DS = Reference('unpublished / DS')
LR = Reference('unpublished / LR')


lineage1         = Lineage('lineage 1',color='magenta', origin='east africa, indian ocean, phillipines')
lineage2         = Lineage('lineage 2',color='blue', origin='TODO')
lineage_beijing  = Lineage('beijing sublineage',lineage2)
lineage3         = Lineage('lineage 3',color='purple', origin='east africa, central asia')
lineage4         = Lineage('lineage 4',color='red', origin='europe, america, africa')
lineage_X        = Lineage('sublineage X',lineage4)
lineage_harlem   = Lineage('sublineage harlem',lineage4)
lineage_ghana    = Lineage('sublineage ghana',lineage4)
lineage_ural     = Lineage('sublineage ural',lineage4)
lineage_vietnam  = Lineage('sublineage vietnam',lineage4)
lineage_LAM      = Lineage('sublineage LAM',lineage4)
lineage_iran     = Lineage('sublineage iran',lineage4)
lineage_uganda   = Lineage('sublineage uganda',lineage4)
lineage_cameroon = Lineage('sublineage cameroon',lineage4)
lineage_PGG3     = Lineage('sublineage PGG3',lineage4)
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

sublineages_SNPs = [

    Test(SNP(genome=ancestor,pos=3798451,base='G'),lineage_X,DS),
    Test(SNP(genome=ancestor,pos= 501615,orig='C',base='G'),lineage_X,DS),
    Test(SNP(genome=ancestor,pos= 918316,base='C'),lineage_X,DS),
    Test(SNP(genome=ancestor,pos=3013784,base='G'),lineage_harlem,DS),
    Test(SNP(genome=ancestor,pos=2759390,base='A'),lineage_harlem,DS),
    Test(SNP(genome=ancestor,pos=3587072,base='C'),lineage_harlem,DS),
    Test(SNP(genome=ancestor,pos=4409231,base='G'),lineage_ghana,DS),
    Test(SNP(genome=ancestor,pos= 108320,base='T'),lineage_ghana,DS),
    Test(SNP(genome=ancestor,pos=3571122,base='A'),lineage_ghana,DS),
    Test(SNP(genome=ancestor,pos=2181026,base='C'),lineage_ural,DS),
    Test(SNP(genome=ancestor,pos=1860528,base='C'),lineage_ural,DS),
    Test(SNP(genome=ancestor,pos=3198496,base='A'),lineage_ural,DS),
    Test(SNP(genome=ancestor,pos=3966059,base='C'),lineage_vietnam,DS),
    Test(SNP(genome=ancestor,pos=1392921,base='A'),lineage_vietnam,DS),
    Test(SNP(genome=ancestor,pos=4238963,base='T'),lineage_vietnam,DS),
    Test(SNP(genome=ancestor,pos=1480024,base='T'),lineage_LAM,DS),
    Test(SNP(genome=ancestor,pos=2097990,base='C'),lineage_LAM,DS),
    Test(SNP(genome=ancestor,pos=1613960,base='T'),lineage_LAM,DS),
    Test(SNP(genome=ancestor,pos=2789341,base='C'),lineage_iran,DS),
    Test(SNP(genome=ancestor,pos=1645599,base='T'),lineage_iran,DS),
    Test(SNP(genome=ancestor,pos=2597696,base='T'),lineage_iran,DS),
    Test(SNP(genome=ancestor,pos= 990626,base='A'),lineage_uganda,DS),
    Test(SNP(genome=ancestor,pos=2440953,base='T'),lineage_uganda,DS),
    Test(SNP(genome=ancestor,pos=1597405,base='A'),lineage_uganda,DS),
    Test(SNP(genome=ancestor,pos=3191099,base='A'),lineage_cameroon,DS),
    Test(SNP(genome=ancestor,pos= 540349,base='T'),lineage_cameroon,DS),
    Test(SNP(genome=ancestor,pos=4260742,base='A'),lineage_cameroon,DS),
    Test(SNP(genome=ancestor,pos=1692141,base='A'),lineage_PGG3,DS),
    Test(SNP(genome=ancestor,pos=1068432,base='A'),lineage_PGG3,DS),
    Test(SNP(genome=ancestor,pos=4255922,base='A'),lineage_PGG3,DS),
]

lineage_11 = Lineage('sublineage 1.1',lineage1)
lineage_12 = Lineage('sublineage 1.2',lineage1)
lineage_13 = Lineage('sublineage 1.3',lineage1)
lineage_13a= Lineage('sublineage 1.3',lineage1, comutant='a')
lineage_13b= Lineage('sublineage 1.3',lineage1, comutant='b')
lineage_14 = Lineage('sublineage 1.4',lineage1)
lineage_14a= Lineage('sublineage 1.4',lineage1, comutant='a')
lineage_15 = Lineage('sublineage 1.5',lineage1)
lineage_15a= Lineage('sublineage 1.5',lineage1, comutant='a')
lineage_16 = Lineage('sublineage 1.6',lineage1)
lineage_16a= Lineage('sublineage 1.6',lineage1, comutant='a')
lineage_16b= Lineage('sublineage 1.6',lineage1, comutant='b')

lineage_31 = Lineage('sublineage 3.1',lineage3)
lineage_32 = Lineage('sublineage 3.2',lineage3)
lineage_33 = Lineage('sublineage 3.3',lineage3)
lineage_34 = Lineage('sublineage 3.4',lineage3)
lineage_34a= Lineage('sublineage 3.4',lineage3, comutant='a')
lineage_34b= Lineage('sublineage 3.4',lineage3, comutant='b')


more_sublineages_SNPs = [

    Test(SNP(genome=ancestor,pos=263057,orig='T',base='G'),lineage_11,LR),
    Test(SNP(genome=ancestor,pos=1944723,orig='G',base='T'),lineage_11,LR),
    Test(SNP(genome=ancestor,pos=4237383,orig='C',base='A'),lineage_11,LR),
    Test(SNP(genome=ancestor,pos=1239798,orig='G',base='T'),lineage_12,LR),
    Test(SNP(genome=ancestor,pos=2679370,orig='G',base='T'),lineage_12,LR),
    Test(SNP(genome=ancestor,pos=4398363,orig='G',base='C'),lineage_12,LR),
    Test(SNP(genome=ancestor,pos=1049161,orig='C',base='A'),lineage_13a,LR),
    Test(SNP(genome=ancestor,pos=2623635,orig='G',base='T'),lineage_13a,LR),
    Test(SNP(genome=ancestor,pos=4326148,orig='C',base='A'),lineage_13a,LR),
    Test(SNP(genome=ancestor,pos=588762,orig='C',base='G'),lineage_13b,LR),
    Test(SNP(genome=ancestor,pos=1264457,orig='A',base='C'),lineage_13b,LR),
    Test(SNP(genome=ancestor,pos=2697916,orig='G',base='C'),lineage_13b,LR),
    Test(SNP(genome=ancestor,pos=1365430,orig='G',base='C'),lineage_13,LR),
    Test(SNP(genome=ancestor,pos=2355363,orig='G',base='T'),lineage_13,LR),
    Test(SNP(genome=ancestor,pos=2751502,orig='G',base='T'),lineage_13,LR),
    Test(SNP(genome=ancestor,pos=2330343,orig='C',base='A'),lineage_13,LR),
    Test(SNP(genome=ancestor,pos=685578,orig='C',base='A'),lineage_14a,LR),
    Test(SNP(genome=ancestor,pos=1353773,orig='C',base='G'),lineage_14a,LR),
    Test(SNP(genome=ancestor,pos=4399494,orig='G',base='C'),lineage_14a,LR),
    Test(SNP(genome=ancestor,pos=15177,orig='C',base='G'),lineage_14,LR),
    Test(SNP(genome=ancestor,pos=877055,orig='A',base='C'),lineage_14,LR),
    Test(SNP(genome=ancestor,pos=4397577,orig='G',base='T'),lineage_14,LR),
    Test(SNP(genome=ancestor,pos=312138,orig='C',base='A'),lineage_15a,LR),
    Test(SNP(genome=ancestor,pos=2091783,orig='C',base='G'),lineage_15a,LR),
    Test(SNP(genome=ancestor,pos=2751072,orig='G',base='T'),lineage_15a,LR),
    Test(SNP(genome=ancestor,pos=313425,orig='G',base='T'),lineage_15,LR),
    Test(SNP(genome=ancestor,pos=2118897,orig='C',base='T'),lineage_15,LR),
    Test(SNP(genome=ancestor,pos=4222852,orig='G',base='C'),lineage_15,LR),
    Test(SNP(genome=ancestor,pos=912464,orig='G',base='C'),lineage_16a,LR),
    Test(SNP(genome=ancestor,pos=1956763,orig='G',base='C'),lineage_16a,LR),
    Test(SNP(genome=ancestor,pos=4134998,orig='G',base='T'),lineage_16a,LR),
    Test(SNP(genome=ancestor,pos=3709403,orig='A',base='C'),lineage_16b,LR),
    Test(SNP(genome=ancestor,pos=1823579,orig='C',base='T'),lineage_16b,LR),
    Test(SNP(genome=ancestor,pos=12234,orig='C',base='G'),lineage_16,LR),
    Test(SNP(genome=ancestor,pos=3070311,orig='C',base='A'),lineage_16,LR),
    Test(SNP(genome=ancestor,pos=2788966,orig='G',base='T'),lineage_16,LR),
    Test(SNP(genome=ancestor,pos=258993,orig='G',base='T'),lineage_31,LR),
    Test(SNP(genome=ancestor,pos=1954269,orig='C',base='A'),lineage_31,LR),
    Test(SNP(genome=ancestor,pos=3504452,orig='G',base='C'),lineage_31,LR),
    Test(SNP(genome=ancestor,pos=2312398,orig='G',base='T'),lineage_32,LR),
    Test(SNP(genome=ancestor,pos=3031358,orig='G',base='C'),lineage_32,LR),
    Test(SNP(genome=ancestor,pos=4161737,orig='C',base='G'),lineage_32,LR),
    Test(SNP(genome=ancestor,pos=1115880,orig='A',base='C'),lineage_33,LR),
    Test(SNP(genome=ancestor,pos=1577705,orig='G',base='T'),lineage_33,LR),
    Test(SNP(genome=ancestor,pos=3175288,orig='G',base='C'),lineage_33,LR),
    Test(SNP(genome=ancestor,pos=1653374,orig='G',base='C'),lineage_34a,LR),
    Test(SNP(genome=ancestor,pos=2306228,orig='C',base='G'),lineage_34a,LR),
    Test(SNP(genome=ancestor,pos=4145150,orig='G',base='C'),lineage_34a,LR),
    Test(SNP(genome=ancestor,pos=586030,orig='C',base='A'),lineage_34b,LR),
    Test(SNP(genome=ancestor,pos=4394539,orig='G',base='C'),lineage_34b,LR),
    Test(SNP(genome=ancestor,pos=1321096,orig='G',base='A'),lineage_34b,LR),
    Test(SNP(genome=ancestor,pos=2306257,orig='T',base='G'),lineage_34,LR),
    Test(SNP(genome=ancestor,pos=3121898,orig='T',base='C'),lineage_34,LR),

]


phylo = PhyloTestsuite(mainlineages_SNPs + beijing_SNPs, VERSION) # + sublineages_SNPs + more_sublineages_SNPs )

