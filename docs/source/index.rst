SIST
====

SIST (Stress-Induced Structural Transitions) is a program for analysing
stress-induced structural transitions in superhelical DNA with a specified
base sequence.

SIST supports calculations for:

* strand separation (melting/SIDD)
* Z-DNA formation
* cruciform extrusion
* competition between melting, Z-DNA, and cruciform transitions

SIST 1.0.0 provides a maintained Conda distribution, an installed ``sist``
command, and regression testing against reference scientific outputs.

Quick start
-----------

Install SIST with Conda:

.. code-block:: bash

   conda install -c ccpbiosim -c conda-forge -c bioconda sist

Run a melting calculation:

.. code-block:: bash

   sist -a M -f sequence.fa

Run a competition calculation and write the result to a file:

.. code-block:: bash

   sist -a A -f sequence.fa -o results.txt

.. toctree::
   :maxdepth: 2
   :caption: Contents

   user-guide
   installation
   source-usage
   development
   citation
