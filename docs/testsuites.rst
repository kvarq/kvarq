
.. _testsuites:

KvarQ testsuites
================

Testsuites define positions to scan for as well as how to interpret mutations.
They have to be loaded (see :ref:`Loading Testsuites <loading-testsuites>` and
the :ref:`GUI settings dialog <settings>`) prior to scanning but also to
analyze ``.json`` data using the :ref:`explorer <explorer>`.  A ``.json`` file
generated with a certain combination of testsuites can only be analyzed using
the explorer if testsuites with the same version are used (this is because
``.json`` files only contain the names of the SNPs but the location within
genes is saved in the testsuites).

KvarQ testsuites can be grouped together in a directory.  These directories can
again be grouped together.  See **for example** the ``testsuites/MTBC``
testsuites that are included in the KvarQ source code (additional testsuites
are linked in as submodules in the ``testsuites/`` directory in the source
distribution or can be `downloaded from github <http://github.com/kvarq>`_):

  - ``testsuites/MTBC/_util.py`` : every file that starts with an underscore
    will **not** be loaded from KvarQ when loading testsuites from a directory,
    but can still be used from other python modules in the same directory via
    ``from _util import ancestor`` (which loads the hypothetical ancestor genome
    from a data file in the same directory
  - ``testsuites/MTBC/phylo.py`` : scans for phylogenetic markers in MTBC
  - ``testsuites/MTBC/resistance.py`` : tests for some common resistance
    mutations in MTBC
  - ``testsuites/MTBC/spoligo.py`` : *in silico* spoligo typing of MTBC
  - ``testsuites/MTBC/legacy/`` : older versions of testsuites mentioned above;
    use these to explore old data files (load same version with which the fastq
    file was originally scanned)


.. _roll-your-testsuite:

Rolling your own testsuite
--------------------------

KvarQ makes it very simple to write new testsuites.  Take as an example the
file included below (can be found in the ``testsuites/`` directory of the
:ref:`source distribution <install-from-source>`).  After having developed
a testsuite and tested it on your data, please `send me a note
<mailto:andreas.steiner@unibas.ch>`_ and I will include a link in the KvarQ
distribution.  See :ref:`ebola14` for a tutorial on how to write a testsuite.

.. literalinclude:: ../testsuites/MTBC/example.py
  :language: python
  :linenos:


.. _testsuites-compatibility:

Versions, Compatibility
-----------------------

The following problems can arise when different versions of testsuites are used

  - A testsuite is not compatible with the KvarQ version that loads it.  To
    avoid this scenario, the module global ``GENES_COMPATIBILITY`` is compared
    with the module global :py:attr:`kvarq.genes.COMPATIBILITY` version of
    KvarQ running it.  The first number must be matched exactly and the second
    number must be equal *or smaller* to the one defined in the genes package.
    Whenever KvarQ introduces new features that break the
    backwards-compatibility with the testsuites, the first number is increased.

  - A testsuite is loaded to display data from a ``.json`` file that was generated
    by testsuite with a different number.  The moduel global ``VERSION``
    tells the version of the testsuite defining it.  Upon a backwards-compatible
    change (e.g. deletion of a previous test), the minor number is increased by one.
    Note that introductions of new tests are **not** backwards compatible because
    the new version of the testsuite will be looking for non-existing tests when
    loading data generated with an old version.

