
#### A [KvarQ] testsuite for ebolavirus

Please see the [KvarQ documentation], especially the [tutorial section] for
more information

All data is taken from:

> "Genomic surveillance elucidates Ebola virus origin and transmission during the
> 2014 outbreak" -- Gire et al (2014) doi 10.1126/science.1259657

The file directory `suppl/` contains files from the [supplementary data]
from that publication, in particular

  - `ebov.mafft.fasta` : sequence alignment of 78 sierra leone strains sequence
    for this paper, 3 previously published guinea strains from 2014, and 
    20 older strains
  - `Table S2 samples.xlsx` name of all samples sequenced in this project;
    used to identify the 78 sierra leone strains
  - `Table S4 2014_specific_snps.xlsx` : contains sheet with all "2014
    specific SNPs"

The file `suppl/_extract_SNPs.py` was used to create `suppl/SNPs.csv` that
lists for every SNP position the number of occurrences of every base for the
different groups (sierraleone, guinea, previous).  `suppl/SNPs.xlsx` is based
on this listing and highlights all the specifici SL1, SL2, and SL3 SNPs from
the paper.  The same script was also used to directly generate the python
code used in the `sierraleone14.py` testsuite.


[KvarQ]: <http://swisstph.ch/kvarq>
[KvarQ documentation]: <http://kvarq.readthedocs.org>
[tutorial section]: <http://kvarq.readthedocs.org/en/latest/tutorial.html>
[supplementary data]: <http://www.sciencemag.org/content/345/6202/1369/suppl/DC1>

