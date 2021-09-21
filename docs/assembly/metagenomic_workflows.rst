====================
Metagenome Assembly
====================

Technical advances in sequencing technologies in recent decades have allowed detailed investigation of complex microbial communities without the need for cultivation, which has proved to be challenging for many communities. Sequencing of microbial DNA extracted directly from enviromental or host-associated samples have provided key information on microbial community composition. These studies have also allowed gene-level characterization of microbiomes, as the first step to understanding the communities functional potential. Furthermore, algorithmic improvements as well as increased availability of computational resources make it now possible to reconstruct whole genomes from metagenomic samples (metagenome-assembled genomes, MAGs). Methods for microbial community composition are discussed in :doc:`../profiling/metagenomes`. Here we describe building :ref:`Metagenomic Assembly`, as well as building :ref:`Gene Catalogs` and :ref:`MAGs` from metagenomic data.

--------------------
Metagenomic Assembly
--------------------

.. image:: ../images/Metagenomic_Assembly_Overview.png


1. **Data Preprocessing**. Before proceeding the the assembly, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`. In addition to standard quality control and adapter trimming, we also suggest normalization with **bbnorm.sh** and merging (See :doc:`../preprocessing/preprocessing` for more details)

2. **Metagenomic Assembly**. Following data preprocessing, we use clean reads to perform a metagenomic assembly using **metaSPAdes**. metaSPAdes is part of SPAdes_ assembly toolkit. assembly. Following the assembly, we generate some assembly statistics using **assembly-stats**, and filter out contigs that are < 1 kbp in length. The script we use for contig filtering can be found here: :download:`contig_filter.py <../scripts/contig_filter.py>`.


.. _SPAdes: https://github.com/ablab/spades

    **Example command**:

    .. code-block:: console

        metaspades.py -t {threads} -m {memory} --only-assembler --pe1-1 <forward_reads.fq.gz> --pe1-2 \
        <reverse_reads.fq.gz> --pe1-s <singeltons.fq.gz> --pe-1m <merged_reads.fq.gz> -o <output_directory>


**Options Explained**

================     =====================================================================================================
-t                   Number of threads.
-m                   Set memory limit in Gb. Spades will terminate if that limit is reached.
--only-assembler     Runs assemlby module only (spades can also perform read error correction, this step will be skipped).
--pe-1               Forward reads.
--pe1-2              Reverse reads.
--pe1-s              Unpaired reads.
--pe-1m              Merged reads.
-o                   Specify output directory.
================     =====================================================================================================

**Example Command for filtering and stats**:

  .. code-block:: console

      python contig_filter.py {params.sample} contigs {sample/contigs.fasta.gz {params.workfolder}/{params.sample}
      assembly-stats -l 500 -t <(zcat {sample.min500.fasta.gz) > {sample}.assembly.stats


.. note::

    **Computational Resources**: Metagenomic assembly requires a lot of memory (usually > 100 Gb).
    You can use multiple threads (16-32) to speed up the assembly


3. The metagenomic contigs generated in step 2 can now be used to build and/or profile :ref:`Gene Catalogs` or to construct :ref:`MAGs`.

--------------
Gene Catalogs
--------------

Gene catalog generation and profiling (i.e. gene abundance estimation) can provide important insights into the community's structure, diversity and functional potential. This analysis could also identify relationships between genetic composition and environmental factors, as well as disease associations.

.. note:: Integrated catologs of reference genes have been generated for many ecosystems (<add link to ocean>, <add link to human gut>), and might be a good starting point for the analysis.


Creation
^^^^^^^^

This protocol will allow you to create a denovo gene catalog from your metagenomic samples.

.. image:: ../images/Building-gene-catalog.png


1. **Gene calling**. We use **prodigal** to extract protein-coding genes from metagenomic assemblies (usually use **contigs** as input). Prodigal has different gene prediction modes with single genome mode as default. To run prodigal on metagenomic mode we add ``-p meta`` option. This will produce a fasta file with amino acid sequences (.faa), nucleotide sequences (.fna), as well as an annotation file (.gff)

    **Example command**:

    .. code-block:: console

        zcat {in.fa.gz} | prodigal -a {out.faa} -d {out.fna} -f gff -o {out.gff} -c -q -p meta

=========    =====================================================================================================
-a
-d
-f
-o
-c
-q
-p
=========    =====================================================================================================


2. **Gene de-replication**. The next step is to remove duplicated sequences from the catalog. (Aggregation across samples?) Called genes are dereplicated using BBTools Dedupe_ and CD-HIT_. Some additional step?

.. _Dedupe: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/

.. _CD-HIT: https://github.com/weizhongli/cdhit/wiki

    **Example command: dereplication**:

    .. code-block:: console

        dedupe.sh -Xmx500G in={in.fasta} out={out.rep.fasta} outd={out.red.fasta} \
        threads=64 absorbrc=f exact=t touppercase=t usejni=t ac=t mergenames=t absorbmatch=t; \

**Options Explained**

=============    =====================================================================================================
-Xmx500G
usejni
in
out
outd
threads
absorbrc
exact
touppercase
ac
mergenames
absorbmatch
=============    =====================================================================================================

    **Example command: clustering**:

    .. code-block:: console

        cd-hit-est -i {out.rep.fasta} -o {out.fasta} -c 0.95 -T 64 \
        -M 0 -G 0 -aS 0.9 -g 1 -r 0 -d 0

=========    =====================================================================================================
-i
-o
-c
-T
-M
-G
-aS
-g
-r
-d
=========    =====================================================================================================


Profiling
^^^^^^^^^

.. image:: ../images/Gene-Catalog-Profiling.png

1. **Read alignment**.
2. **Filtering the alignment files**.
3. **Counting gene abundance**.


-----
MAGs
-----

The Holy Grail of metagenomics is to be able to assemble individual microbial genomes from complex community samples. However assemblies with short read assemblers fails to reconstruct complete genomes. For that reason, binning approaches have been developed to facilitate creation of Metagenome Assembled Genomes (MAGs).


.. image:: ../images/MAGs.png

The first steps (Steps 1 through 3) are the same for MAGs as for :ref:`Gene Catalogs` workflow. This workflow starts with size-filtered metaSPAdes assembled contigs.

1. **All-to-all Alignment**. In this step, quality controlled for each of the metagenomic samples is mapped to each of the metagenomic assemblies using BWA. Map reads from all samples against scaffolds in each other sample. Here we use -a to allow mapping to secondary sites.

    **Example Command**:

    .. code-block:: console

        bwa

.. important::

    **Computational Resources**: !

The generated alignment files are then filtered to only include alignments that are at least 45 nucleotides long, with an identity of >= 97 and covering 80 of the read sequence. The alignment filtering was done using ... Other alternatives?

    **Example Command**:

    .. code-block:: console

        sushicounter

2. **Within- and between-sample abundance correlation for each contig**.

    **Example Command**:

    .. code-block:: console

        metaBAT2

.. note::

    How many samples do I need to benefit?
    Strictly speaking need at least 3, with as few as 20 starting to see improvement in the assemblies

3. **Metagenomic Binning**

    **Example Command**:

    .. code-block:: console

        metaBAT2


4. **Quality Control**. Quality checks: CheckM adn Anvi'o

    Quality Metrics



Taxonomic/Functional annotations -> page for that



Further Reading
^^^^^^^^^^^^^^^
`MetaBat2 Wiki <https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices>`_


Alternative workflow: low abundance metagenome/pooled assembly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
