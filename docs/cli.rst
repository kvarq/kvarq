
.. highlight:: bash

.. _cli:

Using KvarQ Command Line Interface (CLI)
========================================

.. _using-cli:

Using The Command Line
----------------------

Depending on which :ref:`installation instructions <installing>` you
followed, the KvarQ command line utility will be accessible in a
different way

  - installation from source: simply enter KvarQ on the command
    line or alternatively call ``python -m kvarq.cli`` (make sure the
    directory containing ``kvarq/`` is in the ``PYTHONPATH`` if you
    just downloaded & compiled KvarQ without installing it)

  - binary installation windows: go to the directory with the kvarq
    files and start ``kvarq.exe``

  - binary installation OS X: enter the following command in a
    shell: ``/path/to/kvarq.app/Contents/MacOS/python -m kvarq.cli``

In either case you can use all the functionality described below
by using this one comand with appopriate flags.  The different flags
are all described briefly if KvarQ is run with the ``--help`` command
line switch.

Note that most command line arguments apply to a specific **subcommand**, but
some arguments (such as declaration of additional :ref:`testsuites
<testsuites>` directories using the ``-t`` option or logging using the ``-d``
and ``-l`` options) apply to all commands and are therefore specified *before*
the subcommand.


.. _loading-testsuites:

Loading of Testsuites
---------------------

Loading of :ref:`Testsuites <testsuites>` is a two step process:

  1. First all available testsuites are **discovered**.  KvarQ
     looks for testsuites in the following directories 

     a) The directory ``testsuites/`` in the directorythat contains the
        executable.  That is where the testsuites are originally located after
        download.

     b) In the directory ``kvarq_testsuites/`` in the user's home directory.

     c) In the directory ``testsuites/`` in the current working directory

     d) From directories specified with the environment variable
        ``KVARQ_TESTSUITES`` -- use your system's path separator (``;`` on
        windows, ``:`` on most other systems).

     e) From any directories specified with the general ``--testsuite-directory``
        (shorthand ``-t``) command line switch.

     If a testsuite is found several times, the last occurrence is used.  This
     allows for easy modification of existing testsuites: simply copy the files
     into the directory ``./testsuites/`` and modify them.  Because discovery
     takes place later on the current working directory, these modifications
     will override the original testsuites.  This mechanism also easily allows
     to have different versions of the same testsuite -- simply rename the
     directory (e.g. copying ``MTBC/`` to ``MTBC-legacy/`` before applying any
     incompatible changes to the testsuites).

  2. Testsuites are later on loaded from this pool of discovered testsuites
     when necessary.  When :ref:`scanning a .fastq file
     <cli-scan-single-file>`, the testsuites have to be specified explicitly,
     but for most other actions (such as :ref:`showing results
     <cli-illustrate>`), testsuites are loaded automatically.


.. _cli-scan-single-file:

Scanning a File
---------------

The ``scan`` subcommand scans a ``.fastq`` file and saves the results
in a ``.json`` file.  There are many additional parameters to this command.
See the following examples to illustrate some scenarios.

Simplest scenario: Scan a file, showing a progress bar during the scanning
process and save the results (using the whole :ref:`MTBC testsuites
<testsuites-MTBC>`)::

    kvarq scan -l MTBC -p H37v_strain.fastq H37v_strain.json

Being more verbose and copying the log output into a separate file::

    kvarq -d -l kvarq.log scan -l MTBC -p H37v_strain.fastq H37v_strain.json

In the following example, only the phylogenetic testsuite is loaded::

    kvarq scan -l MTBC/phylo H37v_strain.fastq H37v_strain.json

There are many more command line options; in the following example, KvarQ
uses only one thread (this results in a much slower scanning, but would be
advisable if many scans are executed in parallel as scheduled jobs),
specifies explicitly the ``.fastq`` variant (normally this variant is guessed
by peeking into the quality scores), and ignores all reads that are shorter
than 30 base pairs (after quality trimming)::

  kvarq scan -l MTBC -t 1 --variant Sanger -r 30 -p H37v_strain.fastq H37v_strain.json

**During the scanning**, it is possible to obtain some additional statistics by
pressing ``<CTRL-C>`` on the terminal (this does not work on Windows when you
use a MinGW bash prompt). Pressing ``<CTRL-C>`` twice within two seconds will
interrupt the scanning process and proceed to calculate the results with the
data gathered so far.

Usually, the default parameters for quality cut-off and minimum overlap (see
:ref:`configuration-parameters`) work pretty well. If you encounter problems
with a particular ``.fast`` file, refer to the example in
:ref:`determine-scanning-parameters`.


.. _cli-summarize:

Extracting results from a batch of scans
----------------------------------------

Normally, you would run KvarQ over a whole series of ``.fastq`` files
and then in the end extract the relevant information from the resulting
``.json`` files.  The ``summarize`` command allows such an extraction
of summary information from multiple ``.json`` files.  The following
command extracts the results, as reported by the different testsuites,
and saves it to a ``.csv`` file::

  kvarq summarize results/*.json > results.csv


.. _cli-info:

Showing information about testsuites
------------------------------------

The ``info`` commands displays version information and some summary statistics
about testsuites.  Testsuites can be specified the same way as when
:ref:`scanning a file <cli-scan-single-file>`, so this command is handy to
estimate how many templates would be loaded with a given testsuite
selection.  Using the ``-L`` command line switch loads all discovered
testsuites::

  kvarq info -l MTBC
  kvarq info -L


Directly Analysing a .fastq
---------------------------

Use to ``show`` subcommand to analyze ``.fastq`` files directly without
performing a scan of the file.  For example, the readlengths that would result
from a specified quality cutoff can be displayed using::

    kvarq show -Q 13 H37v_scan.fastq


.. _cli-illustrate:

Showing Results
---------------

Some simple analysis of ``.json`` files are possible using the command line,
but the :ref:`GUI explorer <explorer>` is much more powerful.

The subcommand ``illustrate`` can be used to show the final results of
the scanning, as well as detailed information about the coverages or
a histogram of the (quality-cut) readlengths encountered::

    kvarq illustrate -r H37v_strain.json
    kvarq illustrate -c H37v_strain.json
    kvarq illustrate -l H37v_strain.json


.. _cli-update:

Updating Results
----------------

Since a ``.json`` file contains not only the final results but also the
intermediare results (encoded in :py:class:`kvarq.analyse.Coverage`), it is
possible to update the results sections after modifying the code without having
to re-scan the ``.fastq`` file.  The ``.json`` file is updated in-place::

    kvarq -d update H37v_scan.json


.. _cli-more-examples:

More Usage Examples
-------------------

Verify File Format Integrity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check all ``.fastq`` files in a directory structure for file format integrity

.. code-block:: bash

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

