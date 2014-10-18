
.. _testsuites:

KvarQ testsuites
================

Testsuites define positions to scan for as well as how to interpret mutations.
They have to be loaded (see :ref:`Loading Testsuites <loading-testsuites>`)
:ref:`selected prior to scanning <scanner>` but also to analyze ``.json`` data
using the :ref:`explorer <explorer>`.  A ``.json`` file generated with
a certain combination of testsuites can only be analyzed using the explorer if
testsuites with the same version are used (this is because ``.json`` files only
contain the names of the SNPs but the location within genes is saved in the
testsuites).

A KvarQ testsuite is a python source file that defines
a :py:class:`kvarq.genes.Testsuite` with the same name as the python file.
Several of these testsuites can be grouped together within a single directory.
Any number of such directories containing testsuite python files can be stored
in a well defined location from which it is then :ref:`discovered in particular
order <loading-testsuites>`.

For example, the testsuites ``spoligo``, ``resistance``, and ``phylo`` are
grouped together in the directory ``MTBC/`` and can be found in the directory
``testsuites/`` of KvarQ:

  - ``testsuites/MTBC/_util.py`` : every file that starts with an underscore
    will **not** be loaded from KvarQ when loading testsuites from a directory,
    but can still be used from other python modules in the same directory via
    ``from _util import ancestor`` (which loads the hypothetical ancestor genome
    from a data file in the same directory
  - ``testsuites/MTBC/phylo.py`` : scans for phylogenetic markers in MTBC
  - ``testsuites/MTBC/resistance.py`` : tests for some common resistance
    mutations in MTBC
  - ``testsuites/MTBC/spoligo.py`` : *in silico* spoligo typing of MTBC


.. _roll-your-testsuite:

Rolling your own testsuite
--------------------------

KvarQ makes it very simple to write new testsuites.  It is probably easiest to
take a pre-existing testsuite and adapt it to your needs. All testsuites
shipped with KvarQ are well annotated and there are some articles in the
:ref:`tutorial section <tutorial>` that show how to adapt the testsuites
in the ``testsuites/example/`` directory.


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


.. _testsuites-example:

Annotated example
-----------------

The following is a dump of the extensively annotated testsuite
``testsuites/examples/example.py`` included with KvarQ

.. literalinclude:: ../testsuites/examples/example.py
  :language: python

..
  .. _testsuites-MTBC:
  MTBC testsuite
  --------------
  .. _testsuites-examples:
  examples testsuite
  ------------------

