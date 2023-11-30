==================================
Taxonomic Profiling of Metagenomes
==================================

Taxonomic profiling of complex microbial communities is an essential first step in the investigation of relationship between community composition and environmental and/or health factors. The most common approach to community profiling is amplification and classification the 16S rRNA gene. Methods related to 16S rRNA analysis are discussed in detail in :doc:`../profiling/16S`. Recently shotgun metagenomic sequencing has started to replace the amplicon based approaches, as it provides higher resolution information about the microbial community, and resolves some of the biases associated with 16S approach. A number of software tools have been developed to taxonomically profile metagenomic samples. These tools have been benchmarked in `recent studies`_. Here we're going to talk about the use of :ref:`mOTUs` and :ref:`mTAGs` for taxonomic profiling.

.. _recent studies: https://doi.org/10.1101/2021.07.12.451567

--------
mOTUs
--------

    mOTUs_ determines the composition of metagenomic samples using 10 single copy phylogenetic marker genes and an extensive database consisting of reference genomes, metagenomes and metagenome assembled genomes from 23 different environments. Different use cases and applications are discussed in detail in a `recent publication`_ and on `mOTUs website`_. Here we provide a quick reference guide to basic mOTUs functionality.

.. _mOTUs: https://github.com/motu-tool/mOTUs
.. _recent publication: https://doi.org/10.1002/cpz1.218
.. _mOTUs website: https://motu-tool.org/


.. note::

    Please download :download:`sample data <../downloads/profiling.tar.gz>` and :download:`conda environment file <../downloads/profiling.yaml>` for this section if you want to follow along. See the :ref:`tutorials` section for instructions on how to unpack the data and create the conda environment.  mOTUs installation requires database download, so expect it to take a little bit of time.


.. mermaid::

   flowchart LR
        id1( Taxonomic profiling<br/>with mOTUs) --> id2(data preprocessing)
        id2 --> id3(profile<br/>fa:fa-cog mOTUs)
        id3 --> id4(merge<br/>fa:fa-cog mOTUs)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        style id2 fill:#F78A4A,stroke:#F8F7F7,stroke-width:1px
        class id3,id4 tool

1. **Data Preprocessing**. Before taxonomic profiling, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`.


.. important::

    In addition to standard quality control and adapter trimming, we also suggest merging of paired-end reads (see :doc:`../preprocessing/preprocessing` for more details). Using merged reads increases speed and accuracy.


2. **Profile**. Taxonomic profiles for each sample can be generated using mOTUs_ `profile` command. The output profile will consist of identified mOTUs and their abundance.


.. code-block:: bash

    mkdir motus_profiles
    motus profile -f  reads/ERR479298_sub1_R1.fq.gz \
         -r reads/ERR479298_sub1_R2.fq.gz \
         -n ERR479298_sub1 -o motus_profiles/ERR479298_sub1.motus -c -k mOTU -q -p
    motus profile -f  reads/ERR479298_sub2_R1.fq.gz \
         -r reads/ERR479298_sub2_R2.fq.gz \
         -n ERR479298_sub2 -o motus_profiles/ERR479298_sub2.motus -c -k mOTU -q -p


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

    motus merge -i motus_profiles/ERR479298_sub1.motus,motus_profiles/ERR479298_sub1.motus -o motus_profiles/merged.motus

========  ===============================
``-i``     list of mOTU profiles to merge
``-o``     output file name
========  ===============================


--------
mTAGs
--------

mTAGs_ generates taxonomic profiles from short-read metagenomic sequencing data using small subunit of the ribosomal RNA (SSU-rRNA). The mTAGs tool uses a reference database built by clustering sequences within each genus defined in SILVA 138 into OTUs at 97% identity. Each OTU is represented in the database as a degenerate consensus sequence (generated using the IUPAC DNA code). mTAGs_ detects sequencing reads belonging to SSU-rRNA and annotates them through the alignment to consensus reference sequences. For more information about the methods please see the `mTAGs paper`_

.. _mTAGs: https://github.com/SushiLab/mTAGs
.. _mTAGs paper: https://doi.org/10.1093/bioinformatics/btab465


.. mermaid::

   flowchart LR
        id1( Taxonomic profiling<br/>with mTAGs) --> id2(data preprocessing)
        id2 --> id3(profile<br/>fa:fa-cog mTAGs)
        id3 --> id4(merge<br/>fa:fa-cog mTAGs)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        style id2 fill:#F78A4A,stroke:#F8F7F7,stroke-width:1px
        class id3,id4 tool


1. **Data Preprocessing**. As always, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`. As with mOTUs_, we also suggest merging of paired-end reads (see :doc:`../preprocessing/preprocessing` for more details).

2. **Download mTAGs_ database**.

.. code-block:: bash

    mtags download

2. **Profile**. Taxonomic profiles for each sample can be generated using mTAGs_ `profile` command. The tool produces profiles at 8 different taxonomic levels (root, domain, phylum, class, order, family, genus, and otu). Root level combines all domains, the otu level was generated by clustering of sequences within each genus. Each profile will have an 'Unaligned' and 'Unassigned' entry, these represent sequences that could not be aligned or could not be assigned at a given taxonomic level. These need to be taken into account when calculating relative abundances, but should be removed for most of downstream analyses.

.. code-block:: bash

    mkdir mtags_profiles
    mtags profile -f  reads/ERR479298_sub1_R1.fq.gz \
         -r reads/ERR479298_sub1_R2.fq.gz \
         -n ERR479298_sub1 -o mtags_profiles
    mtags profile -f  reads/ERR479298_sub2_R1.fq.gz \
         -r reads/ERR479298_sub2_R2.fq.gz \
         -n ERR479298_sub2 -o mtags_profiles


========  ======================================================================
``-f``      input file(s) for reads in forward orientation
``-r``      input file(s) for reads in reverse orientation
``-s``      input file(s) for unpaired reads (singletons or merged pair end reads)
``-n``      sample name
``-o``      output directory
========  ======================================================================


3. **Merge**. Individual taxonomic profiles can be merged together using  mTAGs_ `merge` on `*.bins` files produced by `mtags profile`.

.. code-block:: bash

    mtags merge -i mtags_profiles/*bins -o mtags_profiles/merged.mtags

========  ===============================
``-i``     list of mOTU profiles to merge
``-o``     output file name
========  ===============================


Choosing between mOTUs and mTAGs
--------------------------------

mOTUs_ and mTAGs_ both generate taxonomic profiles from shotgun metagenomic data, however they differ in their approaches. The choice of the tool will depend on the specific dataset and question at hand.

Here are a few considerations to keep in mind:

#. mTAGs_ and mOTUs_ rely on different methodologies for classification. mTAGs_ uses rRNA sequences clustered at 97% identity, while mOTUs_ relies on 10 universal single-copy marker genes.

#. If you would like to compare your data to rRNA-based studies (for example 16S rRNA amplicon), mTAGs_ would be a better choice.

#. Since mOTUs_ does not rely on rRNA genes (unlike mTAGs_), it avoids the potential problem of copy number variation.

#. mTAGs relies on SILVA database, which in general has a better coverage of diversity. The % of not profiled reads is usually much lower in mTAGs compared to mOTUs. However, this is highly dependent on the environment being studied.

#. Very often the resolution of the mOTUs clusters is higher than that of rRNA OTUs. As a consequence, a single 16S sequence can correspond to multiple mOTUs.

#. The general patterns found in alpha and beta diversity correlate well between these two methods.

#. mOTUs profiles can provide additional information beyond the taxonomic annotation: ref-mOTUs are directly linked to genomes (through specIs defined in ProGenomes2_) and ext-mOTUs are obtained from MAGs. This allows to explore the gene content of the profiled mOTUs, which is not possible for mTAGs profiles, which are defined based on 16S rRNA sequences.

.. _ProGenomes2: https://progenomes.embl.de/

--------
MAPseq
--------

MAPseq_ is a fast and accurate taxonomic classification tool. Since it relies on rRNA sequences for profiling, it can be applied to both amplicon and metagenomic data.

.. _MAPseq: https://doi.org/10.1093/bioinformatics/btx517

.. important::

    Workflow coming soon!
