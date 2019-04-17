
###################################################################
Getting started
###################################################################


Panaurus takes as input pre-computed GFF3 files in the same format as output by
Prokka. These are best produced using the built in script run_prokka.py as this
ensures that the model for calling genes remains consistent between prokka runs.

The script can be invoked by running

.. code-block:: python

   run_prokka -i input_*.fasta -o output_dir --threads 3


The main Panaurus pipeline can then be run as

.. code-block:: python

   panaurus -i output_dir/input_*/*.gff -o output_dir --threads 3


The most important consideration when running Panaurus is how stringent to be
when setting the different thresholds. We provide three default levels via the
'mode' option. One of 'stringent', 'moderate' or 'relaxed'. By default Panaurus
runs in the 'stringent' mode. To run Panaurus with 'relaxed' parameters run


.. code-block:: python

   panaurus -i output_dir/input_*/*.gff -o output_dir --threads 3 --mode relaxed


For a thorough description of all the command line arguments see :doc:`here <command-line-args>`.
