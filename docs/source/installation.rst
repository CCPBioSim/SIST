Installation
============

Supported platforms
--------------------

SIST's Conda package is built and tested on:

* ``linux-64``
* ``osx-64``
* ``osx-arm64``

Conda
-----

The recommended way to install SIST is with Conda:

.. code-block:: bash

   conda install -c ccpbiosim -c conda-forge -c bioconda sist

The package installs the ``sist`` command together with the runtime
dependencies required by SIST, including:

* Python, Biopython, and Beautiful Soup
* Inverted Repeats Finder (IRF)
* the required C++ runtime libraries

SIST 1.0.0 is validated with IRF 3.09, and the Conda package constrains the
runtime dependency to ``>=3.09,<3.10``.

pytest is used for testing and is not required in a normal SIST runtime
environment.

Verify the installation
-----------------------

Confirm that the installed command and runtime dependencies are available:

.. code-block:: bash

   command -v sist
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
* Python 3.12 or later
* IRF 3.09 available as ``irf`` on ``PATH``

Build both C++ components from the repository root:

.. code-block:: bash

   make -C src/trans_three
   make -C src/trans_compete

Install the ``sist`` Python package:

.. code-block:: bash

   python -m pip install .

``sist`` looks for the compiled ``qsidd`` binaries at a Conda install layout
by default, so point it at the binaries just built from source instead:

.. code-block:: bash

   export SIST_TRANS_THREE_BIN="$(pwd)/src/trans_three/qsidd"
   export SIST_TRANS_COMPETE_BIN="$(pwd)/src/trans_compete/qsidd"

The source-tree pipeline can then be run with:

.. code-block:: bash

   sist -a M -f sequence.fa

For cruciform and competition calculations, IRF must be available on ``PATH``.

The source build uses the same calculation modes and command-line parameters as
the installed ``sist`` command. The Source Usage page describes the individual
source components and direct component workflow in more detail.

Deprecated ``master.pl``/``IR_finder.pl`` aliases
--------------------------------------------------

For scripted workflows that still invoke the old Perl script names directly,
``master.pl`` and ``IR_finder.pl`` are also installed as console-script
aliases, backed by the same Python implementation as ``sist``. Both print a
deprecation warning and will be removed in the next release; switch scripts
over to ``sist`` (and the Python API, for direct use of the IR-finder
component) in the meantime.
