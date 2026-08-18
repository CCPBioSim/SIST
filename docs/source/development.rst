Development and Testing
=======================

The maintained development workflow validates both the source tree and the
installed Conda package.

Source build
------------

Build the C++ components from the repository root:

.. code-block:: bash

   make -C trans_three clean
   make -C trans_three
   make -C trans_compete clean
   make -C trans_compete

The Makefiles are also used by the Conda build and accept the compiler and
linker settings supplied by the Conda toolchain.

Source tests
------------

The pytest suite validates command-line behaviour and the maintained scientific
reference outputs.

Python 3.12 or later is supported for source testing. A direct source test run
also requires:

* GNU Make
* a C++ compiler
* Perl
* IRF 3.08 on ``PATH``

Install the Python testing dependencies:

.. code-block:: bash

   python -m pip install -e '.[testing]'

Run the test suite:

.. code-block:: bash

   python -m pytest tests -vv

During a normal source test run, the test fixtures create a temporary copy of
the repository, build the C++ executables, and run the supported calculations
from that working copy.

Scientific regression baselines
--------------------------------

The maintained reference outputs for SIST 1.0.0 are stored under:

.. code-block:: text

   tests/reference/v1.0.0/

The regression suite covers:

* melting
* Z-DNA
* cruciform
* competition

The tests compare deterministic calculation metadata and profile values at the
precision printed by SIST. Runtime is excluded from the scientific comparison.

Reference outputs should only be changed as part of a reviewed scientific
change. A failing regression test should not be resolved by replacing the
reference output without establishing the reason for the difference.

Conda package
-------------

The Conda recipe defines the package build, runtime, and test requirements.

Build requirements
~~~~~~~~~~~~~~~~~~

* Conda C++ compiler
* GNU Make

Runtime requirements
~~~~~~~~~~~~~~~~~~~~

* Perl
* IRF ``>=3.08,<3.09``
* compiler runtime libraries resolved by Conda

Package-test requirements
~~~~~~~~~~~~~~~~~~~~~~~~~

* Python
* pytest

Python and pytest are package-test dependencies only; they are not required for
normal use of the installed SIST package.

Build and test the Conda package
--------------------------------

Build and test the package with:

.. code-block:: bash

   conda build conda-recipe \
       --override-channels \
       -c conda-forge \
       -c bioconda \
       --no-anaconda-upload

``conda-build`` creates isolated build and test environments automatically.

The package test verifies the installed ``sist`` command and its declared
runtime dependencies before running the regression suite.

The two validation paths therefore serve different purposes:

``python -m pytest``
   Builds and tests the maintained source tree.

``conda build``
   Builds the package and tests the installed ``sist`` command using the
   dependencies declared by the Conda recipe.

Release workflow
----------------

The release workflow prepares and publishes a SIST release by:

1. validating the requested version
2. updating the Conda recipe version
3. updating ``CITATION.cff``
4. creating the release tag
5. creating the GitHub release
6. building and testing the Conda package
7. publishing the validated package to the CCPBioSim Anaconda channel

The package published for a release is therefore built and tested from the same
Conda recipe used during development and continuous integration.
