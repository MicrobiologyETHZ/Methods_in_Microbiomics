===================
Data Preprocessing
===================

-----------------------
General Considerations
-----------------------

Data quality control is an essential first step in any bioinformatics workflow. Below we discuss recommended preprocessing steps for **short read Illumina** sequencing data. Broadly, these steps involve Illumina adapter removal, contaminant filtering and quality-trimming. Additional preprocessing steps, recommended only for specific workflows, are detailed in :ref:`Other Considerations`.

.. important::

    This applies to (standard) Illumina short read data. Long read sequencing data from other technologies, or other library preparations from Illumina (e.g. Nextera Mate Pair Reads data) will require a different preprocessing protocol.


.. note::

    Sample data for this section can be found :download:`here <../downloads/Sample1_isolate.tar.gz>`. The conda environment specifications are :download:`here <../downloads/preprocessing.yaml>`. See the :ref:`tutorials` section for intstructions on how to unpack the data and create the conda environment. After unpacking the data, you should have a set of forward (Sample1_R1.fq.gz) and reverse (Sample1_R2.fq.gz) reads. Also included are Illumina adapter sequences (adapters.fa) and PhiX genome (phix174_ill.ref.fa.gz).


.. mermaid::

   flowchart LR
        id1( Preprocessing) --> id2(adapter<br/>trimming<br/>fa:fa-cog BBTools BBDuk)
        id2 --> id3(contaminant<br/>filtering<br/>fa:fa-cog BBTools BBDuk)
        id3 --> id4(quality filtering/<br/>trimming<br/>fa:fa-cog BBTools BBDuk)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        class id2,id3,id4 tool



1.  **Adapter Trimming**. The adapter sequences contain the sequencing primer binding sites, index sequences, and sequences that allow flow-cell binding. Unless removed, these can interfere with downstream analyses. For this and other preprocessing steps, we use  `BBTools <https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/>`_, a set of tools developed by the Joint Genome Institute. Adapter trimming is performed using `BBDuk <https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/>`_. In this step, a FASTA file with Illumina adapter sequences is specified as reference, and BBDuk will perform k-mer matching to trim the adapter sequences from the reads. The example command is shown below.

    **Example command**
        .. code-block:: console

            bbduk.sh -Xmx1G usejni=t in=Sample1_R1.fq.gz in2=Sample1_R2.fq.gz \
            out=Sample1_trimmed_R1.fq.gz out2=Sample1_trimmed_R2.fq.gz \
            outm=Sample1_adapter_matched.fq.gz outs=Sample1_adapter_s.fq.gz  \
            refstats=Sample1.adapter_trim.stats statscolumns=5 overwrite=t ref=adapters.fa \
            ktrim=r k=23 mink=11 hdist=1 2>> preprocessing.log

**Options Explained**

==========    =========================================================================================================
``-Xmx``        This will be passed to Java to set memory usage.
``usejni``      Enable JNI-accelerated version of BBDuk.
``ktrim``       Trims the adapter as well as all the bases to the right of the adapter sequence
``k``           Length of the k-mer used for matching
``mink``        Additionally matches shorter k-mers (with lengths between 23 and 11) to trim partial adapter sequences
``hdist``       Hamming distance for reference k-mers.
``outs``        Write singleton reads whose mate has failed filters to this file.
==========    =========================================================================================================


.. note::

    | `Why are adapter sequences trimmed from only the 3' ends of reads? <https://emea.support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html>`_
    | `Why do we choose k-mer length between 23 and 11? <https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/kmers/kmers>`_

2. **Contaminant removal**. Spike-ins (most commonly PhiX) are usually used for quality control of sequencing runs as well as to ensure nucleotide diversity when sequencing low complexity libraries. We perform this filtering step prior to downstream analysis to be completely sure that these sequences are not be present in your data. Here we use BBDuk and PhiX genome is used as the reference.

    **Example command**

    .. code-block:: console

        bbduk.sh -Xmx1G usejni=t in=Sample1_trimmed_R1.fq.gz in2=Sample1_trimmed_R2.fq.gz \
        out=Sample1_phix_removed_R1.fq.gz out2=Sample1_phix_removed_R2.fq.gz \
        outm=Sample1_phix_matched.fq.gz outs=Sample1_phix_s.fq.gz \
        ref=phix174_ill.ref.fa.gz k=31 hdist=1 \
        refstats=Sample1_phix.stats statscolumns=5 2>> contaminant.log


**Options Explained**

Here, we use a different kmer size

.. note::

    High nucleotide diversity (i.e. equal relative proportions of A,C,G, and T in each cycle) is critical to the performance of Illumina sequencers. Low diversity (or low complexity) libraries, such as amplicon libraries, will have a large proportion of one nucleotide and small proportions of other nucleotides in a cycle. To compensate for low complexity, a PhiX DNA sequence is often added to the library. Different sequencers use different chemistry and image processing software and require different amounts of PhiX spike-in (anywhere between 5% and 50%). Check the latest information about your sequencing platform.


3. **Quality filtering and trimming**. In this step we use BBDuk to trim low quality bases from the ends of the reads and filter reads based on length, average read quality, and number of Ns present.

    **Example command**

    .. code-block:: console

        bbduk.sh -Xmx1G usejni=t in=Sample1_phix_removed_R1.fq.gz in2=Sample1_phix_removed_R2.fq.gz  \
        out1=Sample1_clean_R1.fq.gz out2=Sample1_clean_R2.fq.gz \
        outm=Sample1_qc_failed.fq.gz outs=Sample1_s.fq.gz minlength=45 \
        qtrim=rl maq=20 maxns=1  stats=Sample1_qc.stats statscolumns=5 trimq=14 2>> qc.log

**Options Explained**

================    ==========================================================
``minlength=45``     filters out reads that are shorter than 45 bp
``qtrim=rl``         trims low quality bases on the right and left ends of the reads
``trimq=14``         regions with average quality BELOW 14 will be trimmed
``maq=20``           filters out reads with average quality BELOW 20
``maxns=1``          filters out reads with more than 1 N
================    ==========================================================

.. note::

    Base quality scores (i.e. level of confidence for any one base call) are an integral part of many bioinformatics pipelines (i.e. alignment and variant calling). Quality scores are usually expressed on a Phred scale (:math:`Q=-10log_{10}P`, where P is the probability of an error in the base call). Base quality scores normally range somewhere between 2 and 40, where  Q40 represents an error probability of 1/10000.  More recently, Illumina started using binned quality scores. For example, NovaSeq (with RTA3) only produces 4 Q-scores: 2 is assigned to no-calls, 12 to calls <Q15, 23 to ~Q20 and 37 to >Q30. According to Illumina and in our hands, these binned quality scores did not affect the downstream analyses (i.e. variant calling).


All of the preprocessing commands can be piped together as follows:

.. code-block:: console

    bbduk.sh -Xmx1G usejni=t in=Sample1_R1.fq.gz in2=Sample1_R2.fq.gz \
    out=stdout.fq outm=Sample1_adapter_matched.fq.gz outs=Sample1_adapter_s.fq.gz  \
    refstats=Sample1.adapter_trim.stats statscolumns=5 overwrite=t ref=adapters.fa \
    ktrim=r k=23 mink=11 hdist=1  2>> preprocessing.log | \
    bbduk.sh -Xmx1G usejni=t interleaved=true overwrite=t \
    in=stdin.fq out=stdout.fq outm=Sample1_phix_matched.fq.gz outs=Sample1_phix_s.fq.gz \
    ref=phix174_ill.ref.fa.gz k=31 hdist=1 refstats=Sample1_phix.stats statscolumns=5 2>> preprocessing.log | \
    bbduk.sh -Xmx1G usejni=t overwrite=t interleaved=true \
    in=stdin.fq out1=Sample1_clean_R1.fq.gz out2=Sample1_clean_R2.fq.gz \
    outm=Sample1_qc_failed.fq.gz outs=Sample1_s.fq.gz minlength=45 \
    qtrim=rl maq=20 maxns=1  stats=Sample1_qc.stats statscolumns=5 trimq=14  2>> preprocessing.log;



--------------------
Other Considerations
--------------------

Below are some of the other preprocessing steps that are recommended for specific applications only. All of these steps will be performed on the clean reads produced by general preprocessing workflow outlined above.

========================    ==============================================  ===========
 **Preprocessing Step**               **Recommended for**                    **Tools**
========================    ==============================================  ===========
Filtering out host reads    Any samples containing host DNA                  BBMap
Coverage normalization      Metagenomic assembly (very large samples only)   BBNorm
Paired-read merging         Metagenomic assembly, 16S and mOTUs profiling    BBMerge
========================    ==============================================  ===========

Filtering out host reads
^^^^^^^^^^^^^^^^^^^^^^^^
    Samples containing host DNA can be filtered by mapping the reads to the host genome. This step is performed using `BBMap <https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/>`_ aligner.


.. note::
    Host genome sequences are not provided in the test dataset, but can be downloaded from NCBI, Ensembl, UCSC. Be sure to keep track of the genome version you are using. Genomes for commonly analyzed organisms can also be downloaded from Illumina iGenomes_

.. _iGenomes: https://support.illumina.com/sequencing/sequencing_software/igenome.html

    **Example Command**

    .. code-block::

        bbmap.sh -Xmx23g usejni=t threads=20 overwrite=t qin=33 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast \
        minhits=2 path=host_bbmap_ref qtrim=rl trimq=15 untrim in1=in.1.fq.gz in2=in.2.fq.gz outu1=out.1.fq.gz \
        outu2=out.2.fq.gz outm=out.host.matched.fq.gz 2>> removeHost.log


    This step has to be repeated for singleton sequences generated in the QC step:

    .. code-block::

        bbmap.sh -Xmx23g usejni=t threads=24 overwrite=t qin=33 minid=0.95 maxindel=3 \
        bwr=0.16 bw=12 quickmatch fast    minhits=2 \
        path=host_bbmap_ref qtrim=rl trimq=15 untrim in=in.s.fq.gz outu=out.s.fq.gz \
        outm=out.s.host.matched.fq.gz 2>> out.rmHost.log

=============    ==========================================================
``qin``              Set to 33 or 64 to specify input quality value ASCII offset. 33 is Sanger, 64 is old Solexa. Could be left unspecified (default=auto).
``minid``            Approximate minimum alignment identity to look for.
``maxindel``         Don't look for indels longer than this. Lower is faster.
``bwr``              If above zero, restrict alignment band to this fraction of read length.  Faster but less accurate.
``bw``               Set the bandwidth directly.
``qickmatch``        Generate cigar strings more quickly.
``fast``             Sets other paramters to run faster, at reduced sensitivity.
``minhits``          Minimum number of seed hits required for candidate sites.
``path``             Specify the location to write the index.
``qtrim``            Quality-trim ends before mapping.
``trimq``            Trim regions with average quality below this.
``untrim``           Undo trimming after mapping.
``in``               Primary reads input.
``outu``             Write only unmapped reads to this file.
``outm``             Write only mapped reads, that fail filters to this file.
=============    ==========================================================


Normalization
^^^^^^^^^^^^^
    This step normalizes the coverage by down-sampling reads over high-coverage areas. This step is only necessary for very large metagenomic samples in order to make the assembly computationally tractable. An example using `BBNorm <https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/>`_ is shown below. As above this step needs to be repeated for the singletons.

**Example Command**

    .. code-block::

        bbnorm.sh -Xmx{memory_limit}G threads={threads} extra=s.fq.gz in1=r1.fq.gz \
        in2=r2.fq.gz out1=output_1.fq.gz out2=output_2.fq.gz target=40 mindepth=0 hist=output.hist \
        peaks=output.peaks &> pe_norm.log; \

        bbnorm.sh -Xmx{memory_limit}G threads={threads} extra=r1.fq.gz,r2.fq.gz \
        in=s.fq.gz out=output_s.fq.gz target=40 mindepth=0 hist=output.hist2 \
        peaks=output.peaks2 &> s_norm.log

=============    ==========================================================
``-Xmx``             This will be passed to Java to set memory usage.
``threads``          Set to number of threads desired.
``extra``            For the kmer table: Additional files to use for input, but not for output.
``in1``              Path to the forward reads.
``in2``              Path to the reverse reads.
``out1``             Normalized forward reads.
``out2``             Normalized reverse reads.
``target``           Target normalization depth.
``mindepth``         Kmers with depth below this number will not be included when calculating the depth of a read.
``hist``             Specify a file to write the input kmer depth histogram.
``peaks``            Write the peaks to this file.
=============    ==========================================================

Pair-read Merging
^^^^^^^^^^^^^^^^^

    Merging refers to merging two overlapping reads into one. This is recommended for amplicon data, mOTUs profiling and metagenomic assembly. We do not usually merge the reads for isolate genome assembly. This can be done using `BBMerge <https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/>`_ .

**Example Command**
    .. code-block::

        bbmerge.sh -Xmx32G threads=32 in1=Sample1_R1.fq.gz in2=Sample1_R2.fq.gz out=Sample1.m.fq.gz \
        outu1=Sample1.merge.R1.fq.gz outu2=Sample1.merge.R2.fq.gz minoverlap=16 usejni=t \
        ihist=Sample1.merge.hist &> merge.log

=================     ==========================================================
``-Xmx``               This will be passed to Java to set memory usage.
``threads``            Set to number of threads desired.
``in1``                Path to the forward reads.
``in2``                Path to the reverse reads.
``out``                File for merged reads.
``outu1``              File for forward unmerged reads.
``outu2``              File for reverse unmerged reads.
``minoverlap``         Minimum number of overlapping bases to allow merging.
``ihist``              Insert length histogram output file.
``usejni``             Do overlapping in C code, which is faster.  Requires compiling the C code.
=================     ==========================================================