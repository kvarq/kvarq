
.. _changelog:

Changelog
=========

version 0.12.2
~~~~~~~~~~~~~~

  - new ``examples/`` testsuites and :ref:`tutorials <tutorial>`, leading to
    development of new testsuites ebola-sierraleone14_ and MTBC-SNP-barcodes_
  - new :ref:`testsuite discovering/loading <testsuite-loading>` mechanism
  - easier selection of testsuites via :ref:`GUI <scanner>`
  - better support for :ref:`clonal variants <clonal-variants>`
  - ``MTBC/resistance`` scans ``pncA`` gene; mutations in ``rrsS`` and
    ``rrsK`` (ribosomal RNA) are labeled correctly
  - :ref:`summarize <cli-summarize>` command
  - added support for more ``.fastq.gz`` formats, including paired files
  - increased test coverage

.. _ebola-sierraleone14: https://github.com/kvarq/kvarq-ebola-sierraleone14/archive/master.zip
.. _MTBC-SNP-barcodes: https://github.com/kvarq/kvarq-MTBC-SNP-barcodes/

version 0.12.1
~~~~~~~~~~~~~~

  - support ``.fastq.gz`` files
  - added new ways to :ref:`load testsuites <loading-testsuites>`
  - some bugfixes in ``scripts/table_{scan|combine}.py``
  - moved additional :ref:`testsuites <testsuites>` into their own
    `repositories <https://github.com/kvarq>`_

version 0.11.3
~~~~~~~~~~~~~~

  - (first public version, pushed to github)
  - compile on windows without having to install pthread files
  - moved all testsuites into separate ``testsuite/`` directory
  - introduced :ref:`testsuite compatability checks <testsuites-compatibility>`
  - updated ``MTBC.resistance`` testsuite (added ``katG.279``,
    ``rpoC.N698H``, renamed ``rrsK``)

version 0.11.2
~~~~~~~~~~~~~~

  - renamed to KvarQ (previous name was "pyseq")
  - made xlrd, xlwt optional dependencies (import/export data as .csv)
  - polished GUI somewhat

version 0.11.1
~~~~~~~~~~~~~~

  - internally use :py:class:`kvarq.analyse.Coverage` instead of
    :py:class:`kvarq.genes.Test` to identify sequences and hits
  - more compact file format
  - support legacy ``.json`` files; legacy testsuites can also be loaded
    separately (some are included in ``testsuites/legacy/``)

version 0.10.10
~~~~~~~~~~~~~~~

  - enabled multiple use of same :py:class:`kvarq.genes.Test` for different
    :py:class:`kvarq.genes.Testsuite`

version 0.10.9
~~~~~~~~~~~~~~

  - implemented all plotting in ``Tkinter`` -- ``matplotlib`` not needed anymore
  - added settings dialog to GUI
  - moved non-published tests into separate files in ``testsuites/`` directory

version 0.10.8
~~~~~~~~~~~~~~

  - added new fields to stats : ``nseqhits``, ``records_parsed``
  - added new attributes to :py:class:`kvarq.fastq.Fastq` :
    ``readlength``, ``records_approx``
  - (these fields are also displayed in the json explorer)
  - added terminal color support
  - improved FastQ quality score decoding
  - improved ``scripts/table_combine.py`` (insert data into existing table)
  - include html documentation in distributions
  - testsuites can be loaded from arbitrary locations (see :ref:`roll-your-testsuite`)

version 0.10.7
~~~~~~~~~~~~~~

  - relaxed ``.fastq`` file format specifications

version 0.10.6
~~~~~~~~~~~~~~

  - ``.fastq`` file format checking (both in kvarq.fastq and kvarq.engine)
  - give every thread 10 file junks to keep scanning speed constant until
    the end

version 0.10.5
~~~~~~~~~~~~~~

  - "batch processing" in ``kvarq.gui.{simple|explorer}``
  - resistances output aa number
  - cleaned up output spoligo/resistances
  - do not include hits in .jsons when using gui/simple

