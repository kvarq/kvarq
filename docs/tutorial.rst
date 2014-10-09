
.. _tutorial:

Tutorial
========

The following examples show how KvarQ can be used to quickly analyze genomic
data in ``.fastq`` format.  All examples assume that you have
:ref:`successfully downloaded and installed KvarQ <installing>`, but have no
other prerequisites.  The tutorials are ordered from simpler to more
complicated.


.. _ebola14:

Ebola Outbreak 2014
-------------------

`Gire et al`_ have sequenced some 99 virus genomes during the `2014
Ebola outbreak`_ and immediately released all sequence data to "facilitate
rapid global research".  In the following, we're going to :ref:`develop
a simple testsuite <roll-your-testsuite>` that allows KvarQ to say whether
a ``.fastq`` file is from *a)* a ebola virus, *b)* the 2014 outbreak, and
*c)* from any of the three sublineages defined in the paper (see figure 4A).

.. _Gire et al: http://www.sciencemag.org/content/345/6202/1369.full
.. _2014 Ebola outbreak: https://en.wikipedia.org/wiki/2014_West_Africa_Ebola_virus_outbreak

Creating the Testsuite
~~~~~~~~~~~~~~~~~~~~~~

The following steps led to the `Kvarq Ebola sierraleone14 testsuite`_.

All necessary supplementary materials can be downloaded from the `Science
webpage`_.  In particular, we're interested in

  - Table S2 with the accession numbers of all 99 sequenced genomes.  This
    are serial isolates from 78 patients.

  - File S1 that contains ``ebov.mafft.fasta`` with the sequence alignment.
    There are 20 sequences before 2014 and 81 sequences from 2014.  The
    sequences form Sierra Leone (sequenced in this paper) are identified
    by comparing them to the ID column in Table S2.

  - Table S4 with the sheet ``2014_specific_snps`` that lists all SNPs that
    were found in the new sequences.

With these three files we can generate a list of SNPs that are unique to the
isolates from Sierra Leone.  `Gire et al`_ define three sub-lineages (figure
4A) and by comparing the number of sequences that have specific SNPs we can
compile the following table (see the `script to extract the SNPs`_ on github).

======== ========= ======= ==========
Position Ancestral Derived Sublineage
======== ========= ======= ==========
   800        C        T    SL2
  1849        T        C    SL1
  6283        C        T    SL1
  8928        A        C    SL2
 10218        G        A    SL3
 13856        A        G    SL1
 15660        T        C    SL1
 15963        G        A    SL2
 17142        T        C    SL2
======== ========= ======= ==========

Creating a testsuite from this data is quite straightforward.  First, we choose
a reference genome to extract data from.  This reference genome should have the
ancestral genotype for every SNP that we define.  Because we are only
interested in SNPs from the 2014 strains, we simply take a genome from
a previous isolate, for example the sequence ``EBOV_1976_KC242801`` from the
file ``ebov.mafft.fasta`` and we save it in a new file called ``EBOV76.fasta``.
this file is then loaded as a :py:class:`Genome <kvarq.genes.Genome>` in the
testsuite:

.. code-block:: python

  # old ebola genome from previous outbreak
  EBOV76 = Genome(os.path.join(os.path.dirname(__file__), 'EBOV76.fasta'))

Next, we define the a :py:class:`Reference <kvarq.genes.Reference>` that
identifies the source of the data.  Then, we :py:class:`Genotype
<kvarq.genes.Genotype>` can be bound to a gene, but in our case we simply
specify it by name.

.. code-block:: python

  gire14 = Reference('Gire et al (2014) doi 10.1126/science.1259657')

  # sub-lineages as defined in gire14
  SL1 = Genotype('SL1')
  SL2 = Genotype('SL2')
  SL3 = Genotype('SL3')

In the next step we define the actual :py:class:`SNPs <kvarq.genes.SNP>`
and bind them to the genotypes defined above.

.. code-block:: python

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

Finally, we define a new :py:class:`Testsuite <kvarq.genes.Testsuite>` from
these SNPs and instantiate it to a variable called ``sierraleone14``, which
must be the same name as python file.  We could use the standard testsuite:

.. code-block:: python

  sierraleone14 = Testsuite(SNPs, VERSION)

But we instead choose to define a new testsuite called ``CountGenotype`` that
subclasses the :py:meth:`_analyse <kvarq.genes.Testsuite._analyse>` method, to
summarize all SNPs into one line that shows the genotype and the number of SNPs
found for this genotype.  See `the complete testsuite`_ on github.


Creating a testsuite from this test data is quite straightforward: simply
define each of the :py:class:`SNPs <kvarq.genes.SNP>` as a :py:class:`Test
<kvarq.genes.Test>` and instantiate a :py:class:`Testsuite
<kvarq.genes.Testsuite>`.  As a bonus, the `Kvarq Ebola sierraleone14
testsuite`_ overrides the ``_analyse`` method of the testsuite to display how
many of the specified SNPs have been found for every sublineage.  The reference
genome ``EBOV76.fasta`` is the first genome found in the file
``ebov.mafft.fasta``.

View the complete source code of the testsuite `on github`_.

.. _KvarQ Ebola sierraleone14 testsuite: https://github.com/kvarq/kvarq-ebola-sierraleone14/archive/master.zip
.. _Science webpage: http://www.sciencemag.org/content/345/6202/1369/suppl/DC1
.. _Gire et al: http://www.sciencemag.org/content/345/6202/1369.full
.. _script to extract the SNPs: https://github.com/kvarq/kvarq-ebola-sierraleone14/blob/master/suppl/_extract_SNPs.py
.. _the complete testsuite: https://github.com/kvarq/kvarq-ebola-sierraleone14/blob/master/sierraleone14.py
.. _on github: https://github.com/kvarq/kvarq-ebola-sierraleone14/archive/master.zip

Running the Testsuite
~~~~~~~~~~~~~~~~~~~~~

Choose any of the patients, e.g. ``EM119`` from table S2.  The corresponding
accession number KM233042_ in the nucleotide archive yields another link into
the biosample database : SAMN02951962_, from where the raw sequencing data
can be downloaded.  NCBI stores the ``.fastq`` file in ``.sra`` format, but
this can easily be converted after download using the ``fastq-dump`` command
from the `SRA Toolkit`_. 

Now simply download the `Kvarq Ebola sierraleone14 testsuite`_, start the
:ref:`KvarQ GUI <gui>`, :ref:`load the testsuite <settings>` and analyze the
``.fastq`` file, or launch KvarQ from the :ref:`command line <cli>`.  The
resulting ``.json`` file can be opened in the :ref:`explorer <explorer>` and
should show that the sample from ``EM119`` is sublineage three, showing all the
6 SNPs from SL1, the 4 SNPs from SL2, and the SNP from SL3.

Downloading the sample from ``EM120`` (biosample SAMN02951963_) and analyzing
it the same way shows that this sample also is positive for the 6 SNPs from SL1,
and the 4 SNPs from SL2, but that is missing the SNP from SL3 (opening the
SNP with the explorer shows that it has the original base ``G`` at position
10218).

.. _KM233042: https://www.ncbi.nlm.nih.gov/nuccore/KM233042
.. _SAMN02951962: https://www.ncbi.nlm.nih.gov/biosample/SAMN02951962
.. _SRA Toolkit: http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
.. _SAMN02951963: https://www.ncbi.nlm.nih.gov/biosample/2951963

