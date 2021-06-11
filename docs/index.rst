
.. Microbial Informatics documentation master file, created by
   sphinx-quickstart on Fri Nov 27 14:14:37 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====================
Microbial Informatics
=====================

Set of guidelines and best practices for robust and reproducible bioinformatics processing and data analysis
with the focus on Microbiomics research.


.. important::

    These guidelines are only guidelines. They are constantly under development and subject to change without warning


.. _main-reproducible-research:

---------------------
Reproducible Research
---------------------

Guidelines for creating and maintaining reproducible data workflows


.. _main-bioinformatic-pipelines:

-----------------------
Bioinformatic Pipelines
-----------------------

Collection of bioinformatic pipelines + collection of 'custom scripts' to assist in the analysis


.. _main-data-analysis:

-----------------------
Data Analysis
-----------------------
A forever evolving list of suggestions and considerations when analysing microbiome data

.. _main-support:

-------
Support
-------
* For quesions, comments and suggestions contact ...


---------
Resources
---------


.. toctree::
   :caption: Reproducible Research
   :name: reproducible_research
   :hidden:
   :maxdepth: 1

   reproducible_research/data_management
   reproducible_research/git
   reproducible_research/remote_machines
   reproducible_research/workflow_managers

.. toctree::
   :caption: Bioinformatic Pipelines
   :name: bioinformatic_pipelines
   :hidden:
   :maxdepth: 1

   bioinformatic_pipelines/taxonomic_profiling
   bioinformatic_pipelines/16S
   bioinformatic_pipelines/transcriptomics
   bioinformatic_pipelines/genome_assembly
   bioinformatic_pipelines/functional_profiling
   bioinformatic_pipelines/comparative_genomics
   bioinformatic_pipelines/variant_calling

.. toctree::
   :caption: Data Analysis
   :name: data_analysis
   :hidden:
   :maxdepth: 1

   data_analysis/alpha_diversity
   data_analysis/beta_diversity
   data_analysis/differential_expression
