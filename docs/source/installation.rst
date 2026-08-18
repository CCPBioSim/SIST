Installation
============

Conda
-----

The recommended way to install SIST is with Conda:

.. code-block:: bash

   conda install -c ccpbiosim -c conda-forge -c bioconda sist

The package installs the ``sist`` command together with the runtime
dependencies required by SIST, including:

* Perl
* Inverted Repeats Finder (IRF)
* the required C++ runtime libraries

SIST 1.0.0 is validated with IRF 3.08, and the Conda package constrains the
runtime dependency to ``>=3.08,<3.09``.

Python and pytest are used for testing and are not required in a normal SIST
runtime environment.

Verify the installation
-----------------------

Confirm that the installed command and runtime dependencies are available:

.. code-block:: bash

   command -v sist
   command -v perl
   command -v irf

The commands should resolve inside the active Conda environment.

Running ``sist`` without the required arguments displays the command-line
usage:

.. code-block:: bash

   sist

Building from source
--------------------

SIST can also be built directly from the source repository.

A source build requires:

* a C++ compiler
* GNU Make
* Perl
* IRF 3.08 available as ``irf`` on ``PATH``

Build both C++ components from the repository root:

.. code-block:: bash

   make -C trans_three
   make -C trans_compete

The source-tree pipeline can then be run with:

.. code-block:: bash

   perl master.pl -a M -f sequence.fa

For cruciform and competition calculations, IRF must be available on ``PATH``.

The source build uses the same calculation modes and command-line parameters as
the installed ``sist`` command. The Source Usage page describes the individual
source components and direct component workflow in more detail.
