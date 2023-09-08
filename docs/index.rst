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

.. note::
   For the next section, you have to install conda first. Installation instructions can be found `here <https://docs.conda.io/en/latest/miniconda.html>`_.

For each of the sections we have provided a link to a conda environment file (i.e. a file that specifies which packages to install) and a test dataset that can be used for practice. Download these files and save them to the tutorial directory. Next, you can create and activate the conda environment and extract the test data:

.. code-block::

   conda env create -f <environment name>.yaml
   conda activate <environment name>
   tar -xvzf Sample1_isolate.tar.gz

For  :doc:`/preprocessing/preprocessing`, you need to install `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_. Here are the commands to create the conda environment and unpack the data:

.. code-block::

   mamba env create -f isolate_assembly.yaml
   conda activate isolate_assembly
   tar -xvpf Sample1_isolate.tar.gz

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


---------------
Transcriptomics
---------------
Best practices for transcriptomic and metatranscriptomic data analysis


--------------
Modelling
--------------

Best practices in metagenomic data analysis


-------
Support
-------

* If you have any questions or suggestions leave a comment below!


.. toctree::
   :caption: Data Preprocessing
   :name: preprocessing
   :hidden:
   :maxdepth: 1

   preprocessing/preprocessing
   preprocessing/example12

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
..
   profiling/variants
   profiling/function

..
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
   transcriptomics/metatranscriptomics.rst


.. toctree::
   :caption: Modelling
   :name: modelling
   :hidden:
   :maxdepth: 1

   modelling/modelling.rst

.. toctree::
   :caption: Tn-Seq
   :name: Tn-Seq
   :hidden:
   :maxdepth: 1

   tnseq/tnseq.rst
