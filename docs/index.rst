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

.. _tutorials:

----------
Tutorials
----------

Each of the sections listed below will have a link to a conda environment file and a toy dataset that can be used for practice. After downloading the appropriate files, the conda environment can be created and activated as follows:

.. code-block::

   conda env create -f <example_conda.yaml>
   conda activate example_conda

Data files are provided as a compressed archive and can be extracted with the following command:

.. code-block::

   tar -xvzf <Example.tar.gz>

.. _main-reproducible-research:

---------------------
Data Preprocessing
---------------------

.. warning::

   Before proceeding to any of the bioinformatics workflows, make sure you have good quality data.
   See :doc:`/preprocessing/preprocessing` for more.


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
   :name: preprocessing
   :hidden:
   :maxdepth: 1

   preprocessing/preprocessing

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
   profiling/metagenomes
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

