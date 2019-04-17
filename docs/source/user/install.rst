###################################################################
Installing Panaurus
###################################################################


Installation with brew
========

To do...


Installation with conda
========

To do...


Manual Installation
========

Panaurus relies on a number of third party applications that must be installed
and made available on the command line.

* cd-hit
* Prokka
* Prank
* MAFFT
* Clustal Omega

The pipeline can then be installed by first downloading the code repository and
then running the setup script

.. code-block:: bash
   :linenos:

   git clone https://github.com/gtonkinhill/panaurus
   cd panaurus
   python3 setup.py install
