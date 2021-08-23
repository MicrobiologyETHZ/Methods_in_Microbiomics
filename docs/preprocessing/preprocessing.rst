===================
Data Preprocessing
===================

-----------------------
General Considerations
-----------------------

Data quality control is an essential first step to any bioinformatic workflow. Below we discuss recommended preprocessing steps we recommend for **short read Illumina** sequencing data. Broadly, these steps involve Illumina adapter removal, contaminant filtering and quality-trimming. Additional preprocessing steps, recommended only for specific workflows are detailed in :ref:`Other Considerations`.

.. important::

    This applies to (standard) Illumina short read data. Long read sequencing data from other technologies, or other library preparataions from Illumina (ex. Nextera Mate Pair Reads data) will require a different preprocessing protocol.


.. image:: ../images/Preprocessing.png

1.  **Adapter Trimming**. The adapter sequences contain the sequencing primer binding site, index sequences, and sequences that allow flow-cell binding. Unless removed, these can interfere with downstream analyses. For this and other preproccessing we use Joint Genome Institute developed set of tools `BBTools <https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/>`_. Adpater trimming is perfomed using **BBDuk**. In this step, FASTA file with Illumina adapter sequences is specified as reference, and BBDuk will perform k-mer matching to trim the adapter sequences from the reads. The example command is show below.

    **Example command**
        .. code-block:: console

            bbduk.sh -Xmx1G usejni=t in=<forward_fastq> in2=<reverse_fastq> \
            out=<forward_trimmed_fastq> out2=<reverse_trimmed_fastq> \
            outm=<reads_that_fail_filters> outs=<output_adapter_singletons>  \
            refstats=<output_adapter_stats> statscolumns=5 overwrite=t ref=<input.adapters> \
            ktrim=r k=23 mink=11 hdist=1  2>> <log_file>

**Options Explained**

========    =========================================================================================================
-Xmx1G      sets a limit on memory (?)
usejni=t    makes it faster?
ktrim=r     trim the adapter, as well as all the bases to the right of the adapter sequence
k=23        length of the k-mer used for matching
mink=11     additionally match shorter k-mers (with lengths between 23 and 11)
hdist=1     hamming distance for reference k-mers
========    =========================================================================================================


.. note::

    `Why do we only trim adapters on the 3' ends? <https://emea.support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html>`_

2. **Contaminant removal**. Spike-ins (most commonly PhiX) are usually used for quality control of sequencing runs, as well as to ensure nucleotide diversity when sequencing low complexity libraries. These sequences need to be filtered out prior to downstream analysis.

..note:: PhiX and low complexity libraries: https://emea.support.illumina.com/bulletins/2017/02/how-much-phix-spike-in-is-recommended-when-sequencing-low-divers.html
https://emea.support.illumina.com/bulletins/2016/07/what-is-nucleotide-diversity-and-why-is-it-important.html Check latest information about your sequencing platform.

    Here we also use BBDuk; the command is very similar to the one shown above, however here we provide FASTA file of PhiX genome as the reference.

    **Example command**

    .. code-block:: console

        bbduk.sh -Xmx1G usejni=t in=<forward_fastq> in2=<reverse_fastq> \
        out=<forward_filtered_fastq> out2=<reverse_filtered_fastq> \
        outm={output.phix_matched} outs={output.phix_singletons} \
        ref={phix.fasta} k=31 hdist=1 refstats={output.phix_stats} statscolumns=5 2>> {log.log}


3. **Quality filtering and trimming**. In this step we use BBDuk to trim low quality bases from the ends of the reads, and filter reads based on length, average read quality, and number of Ns present.

    **Example command**

    .. code-block:: console

        bbduk.sh -Xmx1G usejni=t in=<forward_fastq> in2=<reverse_fastq>  \
        fastawrap=10000 out1={output.fq1_clean} out2={output.fq2_clean} \
        outm={output.qc_failed} outs={output.qc_singletons} minlength=45 \
        qtrim=rl maq=20 maxns=1  stats={output.qc_stats} statscolumns=5 trimq=14  2>> {log.log};

**Options Explained**

=============    ==========================================================
minlength=45     filter out reads that are shorter than 45 bp
qtrim=rl         trim low qualit bases on right and left ends of the reads
trimq=14         regions with average quality BELOW 14 will be trimmed
maq=20           filter out reads with average quality below
maxns=1          filter out reads with more than 1 N
=============    ==========================================================

.. note::
    Illumina binned quality scores, doesn't seem to effect downstream analysis pathways in our hands.
    `Illumina binned quality scores <https://www.illumina.com/documents/products/whitepapers/whitepaper_datacompression.pdf>`_. `NovaSeq Quality Scores <https://emea.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf>`_.


All of the preprocessing commands can be piped together as follows:

.. code-block:: console

    bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t in=<forward_fastq> in2=<reverse_fastq> \
    out=stdout.fq outm=<output_adapter_matched> outs=<output_adapter_singletons>  \
    refstats=<output_adapter_stats> statscolumns=5 overwrite=t ref=<input.adapters> \
    ktrim=r k=23 mink=11 hdist=1  2 >> <log_file> | \
    bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f interleaved=true overwrite=t \
    in=stdin.fq out=stdout.fq outm={output.phix_matched} outs={output.phix_singletons} \
    ref={input.phix} k=31 hdist=1 refstats={output.phix_stats} statscolumns=5 2>> {log.log} | \
    bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t overwrite=t interleaved=true \
    in=stdin.fq fastawrap=10000 out1={output.fq1_clean} out2={output.fq2_clean} \
    outm={output.qc_failed} outs={output.qc_singletons} minlength={params.minlen} \
    qtrim=rl maq={params.maq} maxns=1  stats={output.qc_stats} statscolumns=5 trimq={params.trimq}  2>> {log.log};



--------------------
Other Considerations
--------------------

========================    ==============================================  ===========
 **Preprocessing Step**               **Recommended for**                    **Tools**
========================    ==============================================  ===========
Paired-read merging         Metagenomic assembly, 16S and mOTUs profiling
Coverage normalization      Metagenomic assembly
Filtering out host reads    Any samples containing host DNA
========================    ==============================================  ===========


Filtering out host reads
^^^^^^^^^^^^^^^^^^^^^^^^

    **Example Command**
    .. code-block::

        bbmap.sh -Xmx23g usejni=t threads=20 overwrite=t qin=33 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast \
        minhits=2 path={human_bbmap_ref} qtrim=rl trimq=15 untrim in1={in.1.fq.gz} in2={in.2.fq.gz} outu1={out.1.fq.gz} \
        outu2={out.2.fq.gz} outm={out.human.matched.fq.gz} 2>> {out.rmHuman.log}

        # This step has to be repeated for singleton sequences generated in the QC step:

        bbmap.sh -Xmx23g usejni=t threads=24 overwrite=t qin=33 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast    minhits=2 \
        path={human_bbmap_ref} qtrim=rl trimq=15 untrim in={in.s.fq.gz} outu={out.s.fq.gz} outm={out.s.human.matched.fq.gz} 2>> {out.rmHuman.log}

Normalization
^^^^^^^^^^^^^


Pair-read Merging
^^^^^^^^^^^^^^^^^