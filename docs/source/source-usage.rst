Source Usage
============

SIST consists of a Python package, ``sist``, and two C++ implementations of
the transition calculations.

Normal installed use should go through the ``sist`` command. The interfaces
described on this page are useful when building, inspecting, or running the
source tree directly.

Source components
-----------------

``src/sist/``
   The ``sist`` Python package: the command-line interface, orchestration
   that dispatches to the ``qsidd`` binaries, and the Inverted Repeats Finder
   (IRF) integration and energy calculations required for cruciform
   calculations.

``src/trans_three/``
   C++ implementation for analysing strand separation, Z-DNA, and cruciform
   extrusion independently.

``src/trans_compete/``
   C++ implementation for analysing competition between strand separation,
   Z-DNA, and cruciform extrusion.

Running ``sist`` from source
-----------------------------

After building both C++ components and installing the ``sist`` package (see
:doc:`installation`), run the source-tree pipeline with:

.. code-block:: bash

   sist -f <sequence_file> -a <algorithm_type> [options]

For example:

.. code-block:: bash

   sist -a M -f sequence.fa

The available algorithm types are:

* ``-a M``: melting transition only (SIDD)
* ``-a Z``: Z-DNA transition only
* ``-a C``: cruciform transition only
* ``-a A``: competition between melting, Z-DNA, and cruciform transitions

Running ``sist`` without the required arguments displays the available
command-line options.

IRF
---

Cruciform and competition calculations require Inverted Repeats Finder.

``sist``'s IR-finder component invokes ``irf`` from ``PATH``. For source
builds, install a compatible IRF 3.08 executable and ensure that:

.. code-block:: bash

   command -v irf

returns the expected executable.

The maintained Conda package provides this dependency automatically.

Direct C++ usage
----------------

Compile the C++ implementations with:

.. code-block:: bash

   make -C src/trans_three
   make -C src/trans_compete

After compilation, each directory contains a ``qsidd`` executable.

Running ``qsidd`` without the required arguments displays its detailed usage
information. For example, from ``src/trans_three``:

.. code-block:: bash

   ./qsidd -f sequence_file

Cruciform and competition component workflow
--------------------------------------------

When running the components directly, cruciform and competition calculations
require the inverted-repeat energy string that ``sist`` normally computes
internally via :mod:`sist.ir_finder` before invoking ``qsidd -X``. To produce
it directly from Python:

.. code-block:: bash

   python -c "
   from sist.ir_finder import IRFinder
   finder = IRFinder(temperature=310.0, shape='linear')
   print(finder.compute_cruciform_energy_string('sequence_file'))
   "

For a cruciform calculation using ``src/trans_three``:

.. code-block:: bash

   ./qsidd -C -X "string" -f sequence_file

For a competition calculation using ``src/trans_compete``:

.. code-block:: bash

   ./qsidd -X "string" -f sequence_file

Here, ``string`` is the output produced above.

``sist`` coordinates this workflow automatically and is normally the
preferred source-tree entry point.

Working directory
-----------------

For cruciform and competition calculations, :mod:`sist.ir_finder` uses the
basename of the input sequence. The sequence file should therefore be present
in the current working directory when these calculations are run.

Example calculation
-------------------

The repository contains an example competition calculation based on
``pbr322.toy.fa``.

The command is:

.. code-block:: bash

   sist \
       -f pbr322.toy.fa \
       -a A \
       -o pbr322.toy.compete.txt \
       -b \
       -p \
       -r

The example directory also contains IRF intermediate output and an EPS
representation of the competition result.
