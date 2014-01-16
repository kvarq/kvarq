
.. _cli:

Using KvarQ Command Line Interface (CLI)
========================================

.. _using-cli:

Using The Command Line
----------------------

Depending on which :ref:`installation instructions <installing>` you
followed, the kvarq command line utility will be accessible in a
different way

  - installation from source: simply enter ``kvarq`` on the command
    line or alternatively call ``python -m kvarq.cli`` (make sure the
    directory containing ``kvarq/`` is in the ``PYTHONPATH`` if you
    just downloaded & compiled KvarQ without installing it)

  - binary installation windows: go to the directory with the kvarq
    files and start ``kvarq.exe``

..
  does not work (anymore)
  - binary installation OS X: enter the following command in a
    shell: ``/path/to/kvarq.app/Contents/MacOS/python -m kvarq.cli``

In either case you can use all the functionality described below
by using this one comand with appopriate flags.  The different flags
are all described briefly if kvarq is run with the ``--help`` command
line switch.

Note that most command line arguments apply to a specific subcommand, but some
arguments (such as selection of :ref:`testsuites <testsuites>` using the ``-t``
option or logging using the ``-d`` and ``-l`` options) apply to all commands
and are therefore specified **before** the subcommand.


Directly Analysing a .fastq
---------------------------

Use to ``show`` subcommand to analyze ``.fastq`` files directly without
performing a scan of the file.  For example, the readlengths that would result
from a specified quality cutoff can be displayed using::

    kvarq show -Q 13 H37v_scan.fastq


.. _cli-scan-single-file:

Scanning a Single File
----------------------

The ``scan`` subcommand scans a ``.fastq`` file and saves the results
in a ``.json`` file.  There are many additional parameters to this command.
See the following examples to illustrate some scenarios.

Simplest scenario: Scan a file, showing a progress bar during the scanning
process and save the results (using the MTBC :ref:`testsuites <testsuites>`)::

    kvarq -t testsuites/MTBC/ scan -p H37v_strain.fastq H37v_strain.json

Being more verbose and copying the log output into a separate file::

    kvarq -t testsuites/MTBC/ -d -l kvarq.log scan -p H37v_strain.fastq H37v_strain.json

In the following example, only the phylogenetic testsuite is loaded from
the MTBC directory::

    kvarq -t testsuites/MTBC/phylo.py scan H37v_strain.fastq H37v_strain.json

**During the scanning**, it is possible to obtain some additional statistics by
pressing ``<CTRL-C>`` on the terminal (this does not work on Windows when you
use a MinGW bash prompt). Pressing ``<CTRL-C>`` twice within two seconds will
interrupt the scanning process and proceed to calculate the results with the
data gathered so far.

Usually, the default parameters for quality cut-off and minimum overlap (see
:ref:`configuration-parameters`) work pretty well. If you encounter problems
with a particular ``.fast`` file, refer to the examples in
:ref:`determine-scanning-parameters`.


.. _cli-scan-batch-of-files:

Scanning a Batch of Files
-------------------------

The kvarq source distribution comes with some utility scripts that allow to
scan a whole batch of ``.fastq`` files with a single invocation from the
command line. the names of the ``.fastq`` files are taken from a table (which
can be in comma separated values format or alternatively in Microsoft Excel
97/2000/XP/2003 for convenience).  The execution of the script
``scripts/table_scan.py`` will try to find every ``.fastq`` file specified in
the table, then run ``kvarq/cli.py`` and save the resulting ``.json`` files
into the specified output directory::

    python scripts/table_scan.py -f '-t testsuites/MTBC/' selection.xls results/

In the following example, the files in the list are only scanned for
phylogenetic markers and resistance information and the number of threads is
limited to four to save some processing power for the other users while all hit
occurences are recorded for further analysis (the flags specified are exeactly
the same as available for the :ref:`scan subcommand <cli-scan-singlefile>` --
``{logfn}`` will be replaced with a file called ``table_scan.log`` in the
destination directory)::

    python scripts/table_scan.py -f '-l {logfn} -t testsuites/MTBC/ scan -t 4 -H -p' selection.xls results/

In a second step, the original table can be combined with the results of all
the ``.json`` files (the following command will create a new file called
``output/selection.xls`` that contains the original ``selection.xls`` as well
as the result from the scans).  The script can easily be modified to include
other data from the ``.json`` file than the default selection.::

    python scripts/table_combine.py selection.xls results/


.. _cli-illustrate:

Showing Results
---------------

Some simple analysis of ``.json`` files are possible using the command line,
but the :ref:`GUI explorer <explorer>` is much more powerful.

The subcommand ``illustrate`` can be used to show the final results of
the scanning, as well as detailed information about the coverages or
a histogram of the (quality-cut) readlengths encountered::

    kvarq -t testsuites/MTBC illustrate -C H37v_strain.json
    kvarq -t testsuites/MTBC illustrate -r H37v_strain.json
    kvarq -t testsuites/MTBC illustrate -l H37v_strain.json


.. _cli-update:

Updating Results
----------------

Since a ``.json`` file contains not only the final results but also the
intermediare results (encoded in :py:class:`kvarq.analyse.Coverage`), it is
possible to update the results sections after modifying the code without having
to re-scan the ``.fastq`` file.  The ``.json`` file is updated in-place::

    kvarq -t testsuites/MTBC -d update H37v_scan.json


.. _cli-more-examples:

More Usage Examples
-------------------

Verify File Format Integrity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check all ``.fastq`` files in a directory structure for file format integrity::

    #!/bin/bash
    for fastq in `find /tbresearch -name \*.fastq`; do
      python -m kvarq.cli -d show "$fastq" 2>"$0_error.log"
      err="$?"
      echo $err $fastq
      if [ $err -ne 0 ]; then
        # file format error
        base=`basename "$fastq"`
        mv "$0_error.log" "${base%.fastq}.log"
      fi
    done
    rm "$0_error.log"


.. _determine-scanning-parameters:

Determine Scanning Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To find ideal values for the :ref:`configuration-parameters` it's a good idea
to first have a look at the output of ``python kvarq.cli show -Q 13`` (minimum
PHRED score of 13 corresponds to p<0.05).

In the following example, the quality score needs to be lowered to yield
anything useable from the ``.fastq``::

  [   0-   3] 4440 (44%)*****************************************************************
  [   3-   6] 2995 (29%)*******************************************
  [   6-   9] 1221 (12%)*****************
  [   9-  12]  618 ( 6%)*********
  [  12-  15]  364 ( 3%)*****
  [  15-  18]  206 ( 2%)***
  [  18-  21]   93 ( 0%)*
  [  21-  24]   44 ( 0%)
  [  24-  27]   11 ( 0%)
  [  27-  30]    1 ( 0%)
  [  30-  33]    0 ( 0%)
  [  33-  36]    0 ( 0%)
  [  36-  39]    2 ( 0%)
  [  39-  42]    1 ( 0%)
  [  42-  45]    2 ( 0%)
  [  45-  48]    2 ( 0%)

In the next example, the minimum overlap and minimum readlength should be adapted to
something below 25::

  [   0-   2]  183 ( 1%)******
  [   2-   4]  209 ( 2%)*******
  [   4-   6]  611 ( 6%)*********************
  [   6-   8]  839 ( 8%)******************************
  [   8-  10]  896 ( 9%)********************************
  [  10-  12]  822 ( 8%)*****************************
  [  12-  14]  867 ( 8%)*******************************
  [  14-  16]  633 ( 6%)**********************
  [  16-  18]  692 ( 6%)************************
  [  18-  20]  628 ( 6%)**********************
  [  20-  22]  499 ( 5%)*****************
  [  22-  24]  520 ( 5%)******************
  [  24-  26] 1810 (18%)*****************************************************************
  [  26-  28]  706 ( 7%)*************************
  [  28-  30]   82 ( 0%)**

