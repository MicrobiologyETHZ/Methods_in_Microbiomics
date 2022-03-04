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

    This documentation is currently under construction.

.. _tutorials:

----------
Tutorials
----------

Tutorials are available for each section listed below. To start, create a directory for the tutorial of your interest:

.. code-block::

   mkdir <section>_tutorial
   cd <section>_tutorial

For example, for  :doc:`/preprocessing/preprocessing`, run:

.. code-block::

   mkdir preprocessing_tutorial
   cd preprocessing_tutorial

For each of the sections we have provided a link to a conda environment file and a test dataset that can be used for practice. Download these files and save them to the tutorial directory. Next, you can create and activate the conda environment and extract the test data:

.. code-block::

   conda env create -f <environment name>.yaml
   conda activate <environment name>
   tar -xvzf Sample1_isolate.tar.gz

For example, for  :doc:`/preprocessing/preprocessing`, run:

.. code-block::

   conda env create -f preprocessing.yaml
   conda activate preprocessing
   tar -xvzf Sample1_isolate.tar.gz

Now, you are ready to run the example commands given in the corresponding section.

.. _main-reproducible-research:

---------------------
Data Preprocessing
---------------------

.. warning::

   Before proceeding to any of the bioinformatics workflows, make sure you have good quality data.
   See :doc:`/preprocessing/preprocessing` for more.


Download :download:`data preprocessing conda environment file </downloads/preprocessing.yaml>` and :download:`data preprocessing test dataset </downloads/Sample1_isolate.tar.gz>`

-----------------------
Assembly
-----------------------

Best practices for :doc:`/assembly/genome_assembly` and :doc:`/assembly/metagenomic_workflows`.

Download :download:`isolate assembly conda environment file </downloads/isolate_assembly.yaml>` and :download:`isolate assembly test dataset </downloads/Sample1_isolate.tar.gz>`


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

