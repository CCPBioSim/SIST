:html_theme.sidebar_secondary.remove:

SIST
====

.. ccpbiosim-hero::
   :logo-light: _static/logos/SIST-logo-black-text.svg
   :logo-dark: _static/logos/SIST-logo-white-text.svg
   :logo-alt: SIST
   :tagline: Stress-Induced Structural Transitions in superhelical DNA
   :button: user-guide | Get started
   :github: CCPBioSim/SIST

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

.. tab-set::

   .. tab-item:: Install

      Install SIST with Conda:

      .. code-block:: bash

         conda install -c ccpbiosim -c conda-forge -c bioconda sist

   .. tab-item:: Melting

      Run a melting calculation:

      .. code-block:: bash

         sist -a M -f sequence.fa

   .. tab-item:: Competition

      Run a competition calculation and write the result to a file:

      .. code-block:: bash

         sist -a A -f sequence.fa -o results.txt

Explore the documentation
-------------------------

.. grid:: 1 2 3 3
   :gutter: 3

   .. grid-item-card:: 📖 User Guide
      :link: user-guide
      :link-type: doc

      Calculation modes, options and output.

   .. grid-item-card:: 📦 Installation
      :link: installation
      :link-type: doc

      Install SIST with Conda or from source.

   .. grid-item-card:: 🔧 Source Usage
      :link: source-usage
      :link-type: doc

      Run SIST and its parts from a source build.

   .. grid-item-card:: 🧬 Examples
      :link: examples
      :link-type: doc

      Example inputs and outputs to try.

   .. grid-item-card:: 🧪 Development
      :link: development
      :link-type: doc

      Building, testing and releasing SIST.

   .. grid-item-card:: 📝 Citation
      :link: citation
      :link-type: doc

      How to cite SIST in your work.

.. toctree::
   :maxdepth: 2
   :hidden:

   User Guide <user-guide>
   Installation <installation>
   Source Usage <source-usage>
   Examples <examples>
   Development <development>
   Citation <citation>
