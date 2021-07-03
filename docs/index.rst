
.. Microbial Informatics documentation master file, created by
   sphinx-quickstart on Fri Nov 27 14:14:37 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======================
Methods in Microbiomics
=======================

Set of guidelines and best practices for robust and reproducible bioinformatics processing and data analysis
with the focus on Microbiomics research.


.. important::

    This documentation is currently under constructions


.. _main-reproducible-research:

---------------------
Data Preprocessing
---------------------

.. warning::

   Before proceeding to any of the bioinformatics workflows, make sure you have good quality data.
   See :doc:`/data_preprocessing/data_preprocessing` for more.


-----------------------
Assembly
-----------------------

Best practices for assembly of genomes and metagenomes


----------
Profiling
----------

Best practises for profiling of amplicon and metagenomic data


----------
Annotation
----------
SNVs, mVIRs, BGCs.


---------------
Transcriptomics
---------------
Best practices for transcriptomic and metatranscriptomic data profiling


--------------
Modelling
--------------

Best practices in metagenomic data analysis


--------------
Workflows
--------------

Example implementation of best practices


-------
Support
-------
* For quesions, comments and suggestions contact ...


---------
Resources
---------


.. toctree::
   :caption: Data Preprocessing
   :name: data_preprocessing
   :hidden:
   :maxdepth: 1

   data_preprocessing/data_preprocessing

.. toctree::
   :caption: Assembly
   :name: assembly
   :hidden:
   :maxdepth: 1

   assembly/genome_assembly
   assembly/metagenomic_workflows
   assembly/long_read

.. toctree::
   :caption: Profiling
   :name: profiling
   :hidden:
   :maxdepth: 1

   profiling/16S
   profiling/motus
   profiling/mtags
   profiling/variants
   profiling/function

.. toctree::
   :caption: Annotation
   :name: annotation
   :hidden:
   :maxdepth: 1

   annotation/bgc.rst
   annotation/mvirs.rst


.. toctree::
   :caption: Transcriptomics
   :name: transcriptomics
   :hidden:
   :maxdepth: 1

   transcriptomics/transcriptomics.rst


.. toctree::
   :caption: Modelling
   :name: modelling
   :hidden:
   :maxdepth: 1

   modelling/modelling.rst



.. toctree::
   :caption: Workflows
   :name: workflows
   :hidden:
   :maxdepth: 1

   workflows/workflows.rst


.. raw:: html

   <script src="https://utteranc.es/client.js"
           repo="MicrobiologyETHZ/Methods_in_Microbiomics"
           issue-term="pathname"
           theme="github-light"
           crossorigin="anonymous"
           async>
   </script>