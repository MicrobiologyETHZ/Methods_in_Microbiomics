================
Genome Assembly
================

Over the recent years, bacterial whole genome sequencing has become an indispensable tool for microbiologists. While powerful, short read sequencing technologies only allow assembly of draft genomes (i.e. assembly consisting of multiple scaffolds). As illustrated below, during short read assembly individual DNA sequences are first assembled into contigs. The contigs can then be linked together into scaffolds (for example with  `mate-pair reads <https://www.illumina.com/science/technology/next-generation-sequencing/mate-pair-sequencing.html>`_).

Another figure here. Figure text (needs proofreading): DNA is randomly sheared into inserts of known size distribution (i.e., min, max, mean). Paired-end reads may be overlapping providing the possibility to 'merge' reads into one sequence (short-insert paired-end reads (SIPERs). In case they do not overlap (long-insert paired-end reads (LIPERs) or mate pair), the sequence between the paired reads remains unknown, while the length of the inner distance can be estimated due to the known insert size distribution. This information helps the scaffolding process.


.. image:: /images/Genome_Assembly.png


.. note::

    Long read sequencing (with PacBio or Nanopore) offers creation of complete circularized bacterial genomes. However, the bioinformatics methods for this are still in development. They are likely to change as technology develops, and the standard protocols are less well established (See :doc:`/assembly/long_read` for current suggestions).


-----------------------------------------
Isolate genome assembly using short reads
-----------------------------------------

.. image:: /images/Isolate_Genome_Assembly.png

1. **Data Preprocessing**. Before proceeding to the assembly, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`. In addition to standard quality control and adapter trimming, we also suggest normalization with bbnorm.sh and merging (see :doc:`../preprocessing/preprocessing` for more details). Besides the common preprocessing steps, we usually run mOTUs_ on the cleaned sequencing reads, to check for sample contamination or mis-labelling (both occur more frequently than you would expect). For more details please check the :doc:`/profiling/motus` section.

.. _mOTUs: https://github.com/motu-tool/mOTUs

2. **Genome Assembly**. Following data preprocessing, we assemble the cleaned reads using SPAdes_. Subsequently, we generate some assembly statistics using assembly-stats, and filter out contigs that are < 500 bp in length. The script we use for contig filtering can be found here: :download:`contig_filter.py <../scripts/contig_filter.py>`.

.. _SPAdes:
 **Example Command**
.. code-block::

    spades.py -t 4 --isolate --pe1-1 {input.fq1} --pe1-2 {input.fq2} \
    --pe1-s {input.s} -o {params.outdir}

**Options Explained**

================     =====================================================================================================
-t                   Number of threads
--isolate            Use SPAdes isolate mode
--pe-1               Forward reads
--pe1-2              Reverse reads
--pe1-s              Unpaired reads
-o                   Specify output directory
================     =====================================================================================================

**Example Command for filtering and stats**:

  .. code-block:: console

      python contig_filter.py {params.sample} scaffolds {sample/contigs.fasta.gz {params.workfolder}/{params.sample}
      assembly-stats -l 500 -t <(zcat {sample.min500.fasta.gz) > {sample}.assembly.stats



3. **Genome Quality Assesment**. The metrics to evaluate genome quality are calculated using QUAST_. The output will contain information on the number of contigs, the largest contig, total length of the assembly, GC%, N50, L50 and others. If reference genome assembly is available, QUAST_ will also assess misassemblies and try to categorize them.

**N50 and L50**: Given a set of contigs sorted by length in descending order, N50 is the length of the shortest contig added to the lengths of previous contigs to obtain a sum which accounts for at least 50% of the length of the total genome. L50 is the smallest number of contigs that add up to at least 50% of the genome length.

.. image:: /images/n50.png


.. _QUAST: http://quast.sourceforge.net/quast.html

    **Example Command**:

    .. code-block:: console

        quast.py {input.scaf1} -1 {input.r1} -2 {input.r2} -o {params.outDir}



4. Gene Calling and Annotation

 **Example Command**
.. code-block::

    prokka --outdir {params.outdir} --locustag {params.locustag} \
    --compliant --prefix {params.locustag} {input.scaffolds} --force


-----------------------
Alternative Approach
-----------------------
Alternatively





