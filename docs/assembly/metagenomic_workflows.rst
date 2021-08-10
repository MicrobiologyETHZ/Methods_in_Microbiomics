====================
Metagenome Assembly
====================

Technical advances in sequencing technologies in recent decades have allowed detailed investigation of complex microbial communities without the need for cultivation, which has proved to be challenging for many communities. Sequencing of microbial DNA extracted directly from enviromental or host-associated samples have provided key information on microbial community composition. These studies have also allowed gene-level characterization of microbiomes, as the first step to understanding the communities functional potential. Furthermore, algorithmic improvements as well as increased availability of computational resources make it now possible to reconstruct whole genomes from metagenomic samples (metagenome-assembled genomes, MAGs). Methods for microbial community composition are discussed in :doc:`../profiling/metagenomes`. Here we describe building :ref:`Gene Catalogs` and :ref:`MAGs` from metagenomic data.



--------------
Gene Catalogs
--------------

    Gene catalog generation and profiling (i.e. gene abundance estimation) can provide important insights into the community's structure, diversity and functional potential. This analysis could also identify relationships between genetic composition and environmental factors, as well as disease associations.

.. note:: Integrated catologs of reference genes have been generated for many ecosystems (<add link to ocean>, <add link to human gut>), and might be a good starting point for the analysis.


Preprocessing
^^^^^^^^^^^^^

1. Before proceeding the the assembly, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`
2. In addition to standard quality control and adapter trimming, we also suggest normalization with bbnorm.sh and merging (also add link in preprocessing)

Creation
^^^^^^^^

This protocol will allow you to create a denovo gene catalog from your metagenomic samples.

.. image:: ../images/Building-gene-catalog.png


3. Following data preprocessing, we use clean reads to perform a metagenomic assembly using **metaSPAdes**. metaSPAdes is part of SPAdes_ assembly toolkit. assembly. ,

.. _SPAdes: https://github.com/ablab/spades

    **Example command**:

    .. code-block:: console

        metaspades.py -t {threads} -m {memory} --only-assembler --pe1-1 {forward_reads.fq.gz} --pe1-2 \
        {reverse_reads.fq.gz} --pe1-s {singeltons.fq.gz} --pe-1m {merged_reads.fq.gz} -o {spades_output_directory}

.. note::

    **Computational Resources**: Metagenomic assembly requires a lot of memory (usually > 100 Gb).
    You can use multiple threads (16-32) to speed up the assembly



4. Following the assembly, we generate some assembly statistics using **assembly-stats**, and filter out contigs that are < 1 kbp in length.


5. Gene calling.
    We use **prodigal** to extract protein-coding genes from metagenomic contigs. Prodigal has different gene prediction modes with single genome mode as default. To run prodigal on metagenomic mode we add ``-p meta`` option.

    **Example command**:

.. code-block:: console

    metaspades.py -t {threads} -m {memory} --only-assembler --pe1-1 {forward_reads.fq.gz} --pe1-2 \
    {reverse_reads.fq.gz} --pe1-s {singeltons.fq.gz} --pe-1m {merged_reads.fq.gz} -o {spades_output_directory}

6. Gene de-replication

**Example command**:

.. code-block:: console

    metaspades.py -t {threads} -m {memory} --only-assembler --pe1-1 {forward_reads.fq.gz} --pe1-2 \
    {reverse_reads.fq.gz} --pe1-s {singeltons.fq.gz} --pe-1m {merged_reads.fq.gz} -o {spades_output_directory}


Profiling
^^^^^^^^^

.. image:: ../images/Profiling-gene-catalog.png


-----
MAGs
-----

The Holy Grail of metagenomics is to be able to assemble individual microbial genomes from complex community samples. However assemblies with short read assemblers fails to reconstruct complete genomes. For that reason, binning approaches have been developed to facilitate creation of Metagenome Assembled Genomes (MAGs).


The first steps (link to preprocessing) and (link to assemlby) are the same for MAGs as for Gene Catalog workflow. This workflow starts with size-filtered metaSPAdes assembled contigs.

Cross-Mapping
^^^^^^^^^^^^^

**Purpose**:

4. In this step, quality controlled for each of the metagenomic samples is mapped to each of the metagenomic assemblies using BWA. Map reads from all samples against scaffolds in each other sample. Here we use -a to allow mapping to secondary sites.

**Example Command**:

.. code-block:: console
    bwa


5. The generated alignment files are then filtered to only include alingments that are at least 45 nucletides long, with an identity of >= 97 and covering 80 of the read sequence. The alignment filtering was done using ... Other alternatives?

**Example Command**:

.. code-block:: console
    sushicounter

Metagenomic Binning
^^^^^^^^^^^^^^^^^^^

**Purpose:**

6. Metagenomic Binning
    MetaBAT2 for within and between sample coverages for each scaffold
    MetaBAT2 for actual binning

8. Quality checks: CheckM adn Anvi'o

    Quality Metrics

9. Taxonomic/Functional annotations -> page for that

- Strictly speaking need at least 3, with as few as 20 starting to see improvement in the assemblies




Why cross-sample mapping?
^^^^^^^^^^^^^^^^^^^^^^^^

How many samples do I need to benefit?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Further Reading
^^^^^^^^^^^^^^^
`MetaBat2 Wiki <https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices>`_


Alternative workflow: low abundance metagenome/pooled assembly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
