Source Usage
============

SIST consists of a Perl pipeline, an IRF integration script, and two C++
implementations of the transition calculations.

Normal installed use should go through the ``sist`` command. The interfaces
described on this page are useful when building, inspecting, or running the
source tree directly.

Source components
-----------------

``master.pl``
   Pipeline used to run the supported SIST calculations.

``IR_finder.pl``
   Processes Inverted Repeats Finder (IRF) output and produces the inverted
   repeat information required for cruciform calculations, including start
   positions, possible extrusion lengths, and cruciform formation energies.

``trans_three/``
   C++ implementation for analysing strand separation, Z-DNA, and cruciform
   extrusion independently.

``trans_compete/``
   C++ implementation for analysing competition between strand separation,
   Z-DNA, and cruciform extrusion.

Running ``master.pl``
---------------------

After building both C++ components, run the source-tree pipeline with Perl:

.. code-block:: bash

   perl master.pl -f <sequence_file> -a <algorithm_type> [options]

For example:

.. code-block:: bash

   perl master.pl -a M -f sequence.fa

The available algorithm types are:

* ``-a M``: melting transition only (SIDD)
* ``-a Z``: Z-DNA transition only
* ``-a C``: cruciform transition only
* ``-a A``: competition between melting, Z-DNA, and cruciform transitions

Running ``perl master.pl`` without the required arguments displays the
available command-line options.

IRF
---

Cruciform and competition calculations require Inverted Repeats Finder.

``IR_finder.pl`` invokes ``irf`` from ``PATH``. For source builds, install a
compatible IRF 3.08 executable and ensure that:

.. code-block:: bash

   command -v irf

returns the expected executable.

The maintained Conda package provides this dependency automatically.

Direct C++ usage
----------------

Compile the C++ implementations with:

.. code-block:: bash

   make -C trans_three
   make -C trans_compete

After compilation, each directory contains a ``qsidd`` executable.

Running ``qsidd`` without the required arguments displays its detailed usage
information. For example, from ``trans_three``:

.. code-block:: bash

   ./qsidd -f sequence_file

Cruciform and competition component workflow
--------------------------------------------

When running the components directly, cruciform and competition calculations
require the output produced by ``IR_finder.pl``.

Run ``IR_finder.pl`` first:

.. code-block:: bash

   perl IR_finder.pl temperature shape sequence_file

For a cruciform calculation using ``trans_three``:

.. code-block:: bash

   ./qsidd -C -X "string" -f sequence_file

For a competition calculation using ``trans_compete``:

.. code-block:: bash

   ./qsidd -X "string" -f sequence_file

Here, ``string`` is the output produced by ``IR_finder.pl``.

``master.pl`` coordinates this workflow automatically and is normally the
preferred source-tree entry point.

Working directory
-----------------

For cruciform and competition calculations, ``IR_finder.pl`` uses the basename
of the input sequence. The sequence file should therefore be present in the
current working directory when these calculations are run.

Example calculation
-------------------

The repository contains an example competition calculation based on
``pbr322.toy.fa``.

The source-tree command is:

.. code-block:: bash

   perl master.pl \
       -f pbr322.toy.fa \
       -a A \
       -o pbr322.toy.compete.txt \
       -b \
       -p \
       -r

The equivalent installed command is:

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
