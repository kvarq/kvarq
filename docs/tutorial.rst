
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


.. _snp-tutorial:

Creating a new SNP testsuite
----------------------------

After reading the interesting article `A robust SNP barcode for typing
Mycobacterium tuberculosis complex strains
<http://www.nature.com/ncomms/2014/140901/ncomms5812/full/ncomms5812.html>`_
I thought it would be nice to analyze some ``.fastq`` files with that new
barcoding scheme.

To get things done quickly, I was browsing through the testsuites in the
``testsuites/examples`` directory and found a testsuite called ``SNPs.py``
that looked promising.  This testsuite defines a function that loads
SNP declarations from a ``.tsv`` file that can easily be edited with a
popular spreadsheet program.

.. code-block:: python

  here = os.path.dirname(__file__)
  # we borrow the reference from the ../MTBC testsuites
  genome_path = os.path.join(here, os.path.pardir, 'MTBC', 'MTB_ancestor_reference.bases')
  genome = Genome(genome_path, 'MTB ancestor')
  ref = Reference('specify reference here')
  # load SNP information from .tsv file (can be edited with Excel)
  SNPs = tsv2SNPs(os.path.join(here, 'SNPs.tsv'), genome, ref)

The format of the ``.tsv`` file is straightfoward:

========    =======    ===
lineage1    3920109    G/T
lineage1    3597682    C/T
lineage1    1590555    C/T
lineage2    1834177    A/C
lineage2    3304966    G/A
...         ...        ...
========    =======    ===

There is simply an identifier, followed by the position of the SNP within the
reference genome (loaded from the file ``../MTBC/MTB_ancestor_reference.bases``
in the example), then the original base, and finally the derived base.
Actually, the SNPs defined in this example testsuite are the same as the ones
used for the main lineage classification in the testsuite ``MTBC/phylo``.  We
can quickly confirm this by performing a :ref:`scan <cli-scan-single-file>`
using the ``MTBC/phylo`` and the ``examples/SNPs`` testsuites and comparing the
result (type these commands in KvarQ's root directory)

.. code-block:: bash

  kvarq scan -l MTBC/phylo -l examples/SNPs tests/fastqs/N0116_1_hits_1k.fastq N0116_phylo_SNPs.json
  kvarq illustrate -r N0116_phylo_SNPs.json

This should result in the following output::

  examples/SNPs
  -------------
  ['lineage2::SNP1834177AC', 'lineage2::SNP3304966GA']

  MTBC/phylo
  ----------
  'lineage 2 -- low coverage (median below 10x)'

So indeed both testsuites report lineage2 -- because ``examples/SNPs`` does not
subclass :py:class:`kvarq.genes.Testsuite`, the result is simply the list of
SNPs that were found in the file, while ``MTBC/phylo`` fuses the two SNPs into
one lineage result and warns at the same time of low coverage, but that's
material for another tutorial post...

Coming back the SNP barcoding: it's simple enough to compile a list of all SNPs
mentioned in the paper.  It starts like this:

==============    =======     ===
lineage1           615938     G/A
lineage1.1        4404247     G/A
lineage1.1.1      3021283     G/A
lineage1.1.1.1    3216553     G/A
lineage1.1.2      2622402     G/A
lineage1.1.3      1491275     G/A
lineage1.2.1      3479545     C/A
...               ...         ...
==============    =======     ===

So let's first create a new directory for the testsuite-to-be-created, calling
it ``testsuites/MTBC-SNP-barcodes``.  Then we copy the following files

  - ``testsuites/MTBC-SNP-barcodes/coll14.py`` : a copy of the file
    ``testsuites/examples/SNPs.py``, will be modified below

  - ``testsuites/MTBC-SNP-barcodes/coll14.tsv`` : the SNP list extracted from
    the paper; you can `download the list from github`_

.. _download the list from github: https://github.com/kvarq/kvarq-MTBC-SNP-barcodes/blob/master/coll14.tsv


Some parts of the example testsuite have to be modified accordingly

.. code-block:: python
  :emphasize-lines: 1,25,26,27,29
  :linenos:

  VERSION = '0.1'
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
  genome = Genome(os.path.join(here, 'MTB_ancestor_reference_coll.bases'), 'MTB ancestor')
  coll14 = Reference('Coll et al (2014) -- doi: 10.1038/ncomms5812')
  SNPs = tsv2SNPs(os.path.join(here, 'coll14.tsv'), genome, coll14)

  coll14 = Testsuite(SNPs, VERSION)


Remarks

  - line 1 : it doesn't really matter what ``VERSION`` we specify, but it's
    important to increase it when the testsuite is modified to
    :ref:`maintain compatibility <testsuites-compatibility>`

  - line 25 : because the reference genome
    ``MTBC/MTB_ancestor_reference.bases`` that was used in the ``MTBC/phylo``
    testsuite has already the derived base in some of the SNPs defined in
    ``coll14.tsv``, we cannot use it as a reference genome (KvarQ asserts that
    the reference genome has the ancestral base for all defined SNPs to prevent
    errors).  therefore, I have assembled a `new reference genome`_ that has
    the ancestral base for all SNPs

  - line 26 : the reference for the testsuite is the original publication
    from which the SNPs are taken

  - line 27 : the SNPs are read rom the file ``coll14.tsv``

  - line 29 : :ref:`the testsuite must be named like the file
    <testsuites-example>`

.. _new reference genome: https://github.com/kvarq/kvarq-MTBC-SNP-barcodes/blob/master/MTB_ancestor_reference_coll.bases

Now let's see whether KvarQ accepts the new testsuite: the command ``kvarq info
-l MTBC-SNP-barcodes/coll14`` should produce the following output::

  version=0.12.2
  testsuites=MTBC-SNP-barcodes/coll14-0.1[62:3162bp]
  sum=62 tests,3162bp
  sys.prefix=/Library/Frameworks/Python.framework/Versions/2.7

So the new testsuite is accepted and KvarQ tells us that it contains 62 tests
totaling 3162 base pairs (that's 62 times 1 base plus two flanks of 25 base
pairs each).


Running the testsuite
~~~~~~~~~~~~~~~~~~~~~

Let's first scan a single ``.fastq`` to make sure the testsuite works as
expected. For example from the internet: MTB_98_1833_.  Then we scan this file
with our new testsuite::

  kvarq scan -p -l MTBC-SNP-barcodes/coll14 MTB_98_1833.fastq.gz MTB_98_1833.json

.. _MTB_98_1833: ftp:////ftp.broad.mit.edu/pub/annotation/mtuberculosis/diversity/MTB_98_1833.fastq.gz

After a couple of minutes we can examine the result of the scan::

  kvarq illustrate MTB_98_1833.json
  kvarq explorer MTB_98_1833.json

This shows us the file contained at the same time SNPs characteristic for
lineage 2 and lineage 4, and that the reads are quite short (around 35 base
pairs after quality trimming).  Practicaly all SNPs were found (with a coverage
ranging from 20 to 50), most in their ancestral variant.

Ok, so everything seems to work and we can proceed scanning our local library
of ``.fastq`` files, by writing a simple bash script

.. code-block:: bash

  #!/bin/bash
  mkdir results/
  for fastq in /genomes/fastqs/*.fastq; do
      json=`basename "$fastq"`.json
      kvarq -l coll14_scan.log scan -l MTBC-SNP-barcodes/coll14 $fastq results/$json
  done

..
  x*

Some hours later we have scanned for the SNP barcodes of hundreds of genomes,
with a copy of the KvarQ log in the file ``coll14_scan.log`` and a new
``.json`` file in the ``results/`` directory for every genome scanned.  This
information can then be further analyzed using a script that shows all
information in tabular form and can also be `downloaded from github`_ (note
that the script starts with an underscore ``_`` because it is not a testsuite
itself and should not be :ref:`auto-discovered <testsuite-loading>`).

The finished testsuite can also be `found on github`_.

.. _downloaded from github: https://github.com/kvarq/kvarq-MTBC-SNP-barcodes/blob/master/_summarize.py
.. _found on github: https://github.com/kvarq/kvarq-MTBC-SNP-barcodes/

