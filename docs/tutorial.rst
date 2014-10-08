

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

  - File S1 that contains ``ebov.mafft.fasta`` with the sequence alignment
    of 101 strains (of which 81 are from the year 2014)

  - Table S2 with the accession numbers of all 99 sequenced genomes

  - Table S4 with the sheet ``2014_specific_snps``

Looking at all SNPs from table S4 and comparing them with the alignments
yields the following table:

==========  ========  =========  =======
Sublineage  Position  Ancestral  Derived
==========  ========  =========  =======
SL1             2263      C         T
SL1             2314      T         C
SL1             4340      T         C
SL1            10057      A         G
SL1            10065      T         G
SL1            18764      G         A
SL2              800      C         T
SL2             8928      A         C
SL2            15963      G         A
SL2            17142      T         C
SL3            10218      G         A
==========  ========  =========  =======


Creating a testsuite from this test data is quite straightforward: simply
define each of the :py:class:`SNPs <kvarq.genes.SNP>` as a :py:class:`Test
<kvarq.genes.Test>` and instantiate a :py:class:`Testsuite
<kvarq.genes.Testsuite>`.  As a bonus, the `Kvarq Ebola sierraleone14
testsuite`_ overrides the ``_analyse`` method of the testsuite to display how
many of the specified SNPs have been found for every sublineage.  The reference
genome ``EBOV76.fasta`` is the first genome found in the file
``ebov.mafft.fasta``.

View the complete source code of the testsuite `on github`_.

.. _on github: https://github.com/kvarq/kvarq-ebola-sierraleone14/blob/master/sierraleone14.py
.. _KvarQ Ebola sierraleone14 testsuite: https://github.com/kvarq/kvarq-ebola-sierraleone14/archive/master.zip
.. _Science webpage: http://www.sciencemag.org/content/345/6202/1369/suppl/DC1

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

