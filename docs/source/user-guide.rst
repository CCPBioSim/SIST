User Guide
==========

Purpose of SIST
---------------

SIST analyses three types of structural transition in superhelical DNA
molecules of specified base sequence and length:

* strand separation
* B-Z transitions
* cruciform extrusion

The statistical-mechanical methods and algorithms used in these analyses are
described in the publications listed on the Citing SIST page.

Running SIST
------------

When installed through Conda, SIST is run with the ``sist`` command.

The general form is:

.. code-block:: bash

   sist -f <sequence_file> -a <algorithm_type> [options]

Input sequences
---------------

A sequence file is supplied with the ``-f`` option. FASTA format is
recommended.

A single input file must contain only one sequence. SIST converts the supplied
sequence into the format required by the selected calculation.

For cruciform and competition calculations, the input sequence should be in
the current working directory and supplied by file name rather than by a path.
This is required by the IRF integration used by SIST 1.0.0.

Calculation modes
-----------------

``-a M``
   Melting transition only (SIDD).

``-a Z``
   Z-DNA transition only.

``-a C``
   Cruciform transition only.

``-a A``
   Competition between melting, Z-DNA, and cruciform transitions.

Cruciform and competition calculations use Inverted Repeats Finder (IRF).
When SIST is installed through Conda, IRF is installed automatically as a
runtime dependency.

Examples
--------

Melting
~~~~~~~

Run a melting calculation using the default parameters:

.. code-block:: bash

   sist -a M -f sequence.fa

Z-DNA
~~~~~

Run a Z-DNA calculation:

.. code-block:: bash

   sist -a Z -f sequence.fa

Cruciform
~~~~~~~~~

Run a cruciform calculation:

.. code-block:: bash

   sist -a C -f sequence.fa

Competition
~~~~~~~~~~~

Run a competition calculation using the default parameters:

.. code-block:: bash

   sist -a A -f sequence.fa

Write the output to a file:

.. code-block:: bash

   sist -a A -f sequence.fa -o output.txt

Run the competition calculation for a circular plasmid at a superhelical
density corresponding to -0.07 and include the additional base, parameter, and
ensemble-average output:

.. code-block:: bash

   sist -a A -s 0.07 -c -b -p -r -f sequence.fa

Command-line options
--------------------

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Option
     - Description
   * - ``-f FILE``
     - Required. Specify the input sequence file.
   * - ``-a MODE``
     - Required. Select ``M``, ``Z``, ``C``, or ``A``.
   * - ``-T VALUE``
     - Set the temperature in kelvin. The default is 310 K.
   * - ``-s VALUE``
     - Set the superhelical density. The default calculation corresponds to
       -0.06. The command-line convention uses a value such as ``0.07`` for a
       reported stress level of -0.07.
   * - ``-i VALUE``
     - Set the ionic strength. The default is 0.01 M.
   * - ``-th VALUE``
     - Set the energy threshold. The default is 12 kcal/mol.
   * - ``-c``
     - Treat the molecule as circular. The default is linear.
   * - ``-n``
     - Use nearest-neighbour melting energetics. The default is copolymeric
       melting energetics.
   * - ``-b``
     - Include the base at each sequence position in the output.
   * - ``-p``
     - Include calculation parameters in the output.
   * - ``-r``
     - Include ensemble-average results in the output.
   * - ``-o FILE``
     - Write the selected SIST output to ``FILE``.

.. note::

   The SIST 1.0.0 parser accepts uppercase ``-T`` for temperature.

Output
------

Without ``-o``, the selected SIST output is written to standard output.

With ``-o FILE``, the selected output is written to the specified file:

.. code-block:: bash

   sist -a M -f sequence.fa -o melting.txt

The exact output sections depend on the calculation mode and selected options.
The ``-b``, ``-p``, and ``-r`` options add base-pair, parameter, and
ensemble-average information respectively.

IRF may print progress information while cruciform or competition calculations
are running.

Algorithm limitations and recommended parameter ranges
------------------------------------------------------

SIST reports warnings when calculations are performed outside the recommended
parameter ranges.

.. warning::

   Calculations outside these ranges may be inaccurate, may fall outside
   physiological conditions, or may require substantially longer execution
   times.

1. The input sequence length should be greater than 1,500 base pairs and less
   than 10,000 base pairs.

2. Energy thresholds (``-th``) below 9 may yield inaccurate results, while
   thresholds above 15 may result in very long execution times.

3. An absolute superhelical density greater than 0.15 may be outside the
   physiological range.

4. Temperatures below 220 K or above 320 K may be outside the physiological
   range.

5. Salt concentrations below 0.0001 M may be outside the physiological range.
