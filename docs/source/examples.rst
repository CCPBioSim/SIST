Examples
========

Competition calculation
------------------------

Example of the output generated when the competition algorithm analyses the
``pbr322.toy.fa`` sequence. These results are described in the Bioinformatics
citation (see :doc:`citation`).

Input file
~~~~~~~~~~

:download:`pbr322.toy.fa <_static/examples/pbr322.toy.fa>`

Command
~~~~~~~

.. code-block:: bash

   sist -f pbr322.toy.fa -a A -o pbr322.toy.compete.txt -b -p -r

Output files
~~~~~~~~~~~~

* :download:`one_line.pbr322.toy.fa <_static/examples/one_line.pbr322.toy.fa>`,
  the sequence file format required by Inverted Repeat Finder (IRF)
* :download:`one_line.pbr322.toy.fa.2.10.10.80.10.20.10000.100.1.html <_static/examples/one_line.pbr322.toy.fa.2.10.10.80.10.20.10000.100.1.html>`,
  an IRF output file
* :download:`one_line.pbr322.toy.fa.2.10.10.80.10.20.10000.100.1.txt.html <_static/examples/one_line.pbr322.toy.fa.2.10.10.80.10.20.10000.100.1.txt.html>`,
  an IRF output file
* :download:`pbr322.toy.compete.txt <_static/examples/pbr322.toy.compete.txt>`,
  the SIST output file

Figure
~~~~~~

:download:`all.eps <_static/examples/all.eps>` represents graphically the
output shown in ``pbr322.toy.compete.txt``.
