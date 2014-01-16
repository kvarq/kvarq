
.. _hacking:

Hacking KvarQ
=============

This chapter is a good starting point if you intend to change anything in
the KvarQ source code, especially in the modules :py:mod:`kvarq.genes` and
:py:mod:`kvarq.analyse`. If you're merely interested in writing a new
:ref:`testuite <testsuite>`, you might also find some valuable information
here in addition to the example testsuites shipped with KvarQ.  Make sure
to have read the :ref:`overview <overview>` first.


.. _about-flanks:

About flanks
------------

In general, KvarQ searches for a given sequence within the genome and checks
for mutations within that "region of interest".  Because ``.fastq`` reads are
only accepted if the overlap between the target sequence and the read has a
specified minimum overlap (:ref:`minreadlength <configuration-parameters>`),
the coverage rapidly falls on both extremities of the target sequence.  For
this reason, a region of interest is flanked with supplementary bases on either
side.  Mutations in these flanks are disregarded, as their sole purpose is
to avoid to have the border effect with low coverage within the region of
interest.  :py:class:`TemplateFromGenome <kvarq.genes.TemplateFromGenome>`
read their bases from a :py:class:`Genome <kvarq.genes.Genome>` and can easily
add an arbitrary flank on either side of the region of interest.


.. _about-positions:

About positions
---------------

  - All positions within a :py:class:`Genome <kvarq.genes.Genome>` are relative
    to the index **one**.  I.e.  the first codon of the genome (if it were
    coding) would be the bases 1-3 which correspond to the first three bytes in
    the file (at file positions 0-2).  Positions increase in the reading direction
    of the plus strand.

  - The :py:class:`TemplateFromGenome <kvarq.genes.TemplateFromGenome>`'s
    ``start`` and ``stop`` attributes refer to positions in the genome.

  - :py:class:`Sequences <kvarq.genes.Sequence>` are simply strings of bases
    and therefore, indexing within sequences starts at the first character of
    this sequence string.  This sequence string may contain flanks on both
    sides.  The (optional) attribute :py:attr:`start
    <kvarq.genes.Sequence.start>` refers to the first base **after** the flank
    (i.e. the first base of the sequence of interest).

  - In a :py:class:`Coverage <kvarq.analyse.Coverage>`, positions simply refer
    to its :py:class:`Sequences <kvarq.genes.Sequence>` on the plus strand. A
    coverage has no ``pos`` attribute and the ``start``, ``stop`` attributes
    correspond to the (plus strand) sequence's ``left`` and ``right``.

  - The :py:class:`Hit <kvarq.engine.Hit>` structure


.. _about-strands:

About the complementary DNA strand
----------------------------------

In general, everything refers to the ``+`` strand of the genome.  Just before
scanning, the analyser creates a complementary copy of every
:py:class:`Sequence <kvarq.genes.Sequence>` in :py:meth:`Analyser.scan()
<kvarq.analyse.Analyser.scan>`.  When assembling the hits in
:py:meth:`Coverage.apply_hit() <kvarq.analyse.Coverage.apply_hit>`, hits on
the complementary strand are mapped on the positive strand again using the
methods :py:meth:`Sequence.plus_idx() <kvarq.genes.Sequence.plus_idx>` and
:py:meth:`Sequence.plus_base() <kvarq.genes.Sequence.plus_base>`.  Finally,
the :py:class:`TemplateFromGenome <kvarq.genes.TemplateFromGenome>` and the
:py:class:`Gene <kvarq.genes.Gene>` attached to a :py:class:`Test
<kvarq.genes.Test>` know whether the ``+`` or the complementary ``-`` strand
is coding and accordingly converts the base mutations in (non-) synonymous
amino acid changes.


.. _sequence-of-tests:

Sequence of tests
-----------------

The :py:class:`Analyser <kvarq.analyser.Analyser>` is initialized with a set of
:ref:`testsuites <testsuites>` that define each a given number of
:py:class:`Test <kvarq.genes.Test>`.  Just before scanning, the analyser
creates a list (actually a :py:class:`OrderedDict <collections.OrderedDict>`)
that contains every test of every testsuite.  The base sequences the
:py:mod:`engine <kvarq.engine>` searches in the ``.fastq`` version is ordered
in the same sequence (with the complementary base sequences added at the end of
the list).  Later on, the :py:class:`Coverage <kvarq.analyse.Coverage>` are
ordered in the same sequence and everything that gets saved to the ``.json``
file (in :py:meth:`encode() <kvarq.analyse.Analyser.encode>`) uses the same
sequence again.  To be able to reconstruct the same order when loading a
``.json`` file, the :py:meth:`decode() <kvarq.analyse.Analyser.decode>` tries
to identify every test by its name (as returned by its :py:meth:`__str__()
<kvarq.genes.Test.__str__>` method) and reconstructs the same sequence of
tests.  The :py:meth:`Analyser[...] <kvarq.analyse.Analyser.__getitem__>`
method retunrs a :py:class:`Coverage <kvarq.analyse.Coverage>` and accepts
integers as well as strings and :py:class:`Test <kvarq.genes.Test>`).


.. _clonal-variants:

About clonal variants ("heterozygous calls")
--------------------------------------------

  - phylogenetic : in current implementation >50% must be present of SNP to
    be accepted; this leads to complete dominance of the clone with >50%
    genetic material in the mixture; see :py:meth:`kvarq.genes.SNP.validate`

