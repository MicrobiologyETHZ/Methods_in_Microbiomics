==================================
Taxonomic Profiling of Metagenomes
==================================

Taxonomic profiling of complex microbial communities is an essential first step in the  investigation of relationship between community composition and environmental and/or health factors. The most common approach to community profiling is amplification and classificaton the 16S rRNA gene. Methods related to 16S rRNA analysis are discussed in detail in :doc:`../profiling/16S`. Recently shotgun metagenomic sequencing has started to replace the amplicon based approaches, as it provides higher resolution information about the microbial community, and resolves some of the biases associated with 16S approach. A number of software tools have been developed to taxonomically profile metagenomic samples. These tools have been profiled in detail in recent studies<add link>. Here we're going to talk about the use of :ref:`mOTUs` and :ref:`mTAGs` for taxonomic profiling.

--------
mOTUs
--------

    mOTUs_ determines the composition of metagenomic samples using 10 single copy phylogenetic marker genes and an extensive database consisting of reference genomes, metagenomes and metagenome assembled genomes from 23 different environments. Different use cases and applications are discussed in detail in our `recent publication`_ and on `mOTUs website`_. Here we provide a quick reference guide to basic mOTUs functionality.

.. _mOTUs: https://github.com/motu-tool/mOTUs
.. _recent publication: https://doi.org/10.1002/cpz1.218
.. _mOTUs website: https://motu-tool.org/


.. note::

    Sample data and conda environment file for this section can be found :download:`here <../downloads/metag.tar.gz>`. See the :ref:`tutorials` section for instructions on how to unpack the data and create the conda environment. ... Installation requires database download expect it to take a little bit of time.


.. mermaid::

   flowchart LR
        id1( Taxonomic profiling<br/>with mOTUs) --> id2(data preprocessing)
        id2 --> id3(profile<br/>fa:fa-cog mOTUs)
        id3 --> id4(merge<br/>fa:fa-cog mOTUs)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        style id2 fill:#F78A4A,stroke:#F8F7F7,stroke-width:1px
        class id3,id4 tool

|
1. **Data Preprocessing**. Before proceeding to the assembly, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`. In addition to standard quality control and adapter trimming, we also suggest merging of paired-end reads (see :doc:`../preprocessing/preprocessing` for more details). Using merged reads increases speed and accuracy.


2. **Profile**. Taxonomic profiles for each sample can be generated using mOTUs_ `profile` command. The output profile will consist of identified mOTUs and their abundance.


.. code-block:: bash

    motus profile -f ERR479298_sample.1.fq.gz  -r ERR479298_sample.2.fq.gz \
    -s ERR479298_sample.s.fq.gz,ERR479298_sample.m.fq.gz -n {sample} \
    -o motus_profiles -c -k mOTU -q -p


========  ======================================================================
``-f``      input file(s) for reads in forward orientation
``-r``      input file(s) for reads in reverse orientation
``-s``      input file(s) for unpaired reads (singletons or merged pair end reads)
``-n``      sample name
``-o``      output file name
``-c``      print result as counts instead of relative abundances
``-k``      taxonomic level (kingdom, phylum, class, order, family, genus, mOTU)
``-q``      print the full rank taxonomy
``-p``      print NCBI taxonomy identifiers
========  ======================================================================


.. important::
    Expect mOTU counts (when run with ``-c`` option) to be relatively small (compared to total number of reads in your sample). The counts are proportional to the library size, and you can expect ~600 mOTU counts for 5,000,000 reads. If you still think you should be getting higher counts, please see FAQ_ for common issues.

.. _FAQ: https://github.com/motu-tool/mOTUs/wiki/FAQ

.. note::

    The unassigned at the end of the profile file represents the fraction of unmapped reads. This represents species that we know to be present in the sample, but we are not able to quantify individually; hence we group them together into an unassigned fraction. For almost all the analysis, it is better to remove this value, since it does not represent a single species/clade. Please see FAQ_ for more information.


3. **Merge**. Individual taxonomic profiles can be merged together using  mOTUs_ `merge` command to facilitate downstream analysis.

.. code-block:: bash

    motus merge -i ERR479298_sample.motus,ERR479298_sample.motus -o merged.motus \
    -a bioreactor,bee,freshwater,human,marine,mouse,soil,wastewater

========  ===============================
``-i``     list of mOTU profiles to merge
``-o``     output file name
``-a``     add pre-computed profiles from different environmental samples: air, bioreactor, bee, cat, cattle, chicken, dog, fish, freshwater, human, marine, mouse, pig, sheep, soil, termite, wastewater
========  ===============================



--------
mTAGs
--------

.. note::
    Differences between mTAGs and mOTUs.

