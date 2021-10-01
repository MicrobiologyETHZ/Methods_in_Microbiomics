================
Genome Assembly
================

Over the recent years, bacterial whole genome sequencing has become indispensible tool for microbiologists. While powerful, short read sequencing technologies only allow assembly of draft genomes (i.e. ). On the other hand, long read sequencing offers creation of complete circularized genomes, however the bioinformatic methods for this are still in development, and are likely to change as technology develops, and the standard protocols are less well established.


-----------------------------------------
Isolate genome assembly using short reads
-----------------------------------------

1. Data Preprocessing.
2. Assembly

 **Example Command**
.. code-block::

    spades.py -t 4 --isolate --pe1-1 {input.fq1} --pe1-2 {input.fq2} \
    --pe1-s {input.s} {params.merged} -o {params.outdir}

3. (Optional) Check for contamination. Link to mOTUs page.


4. Gene Calling and Annotation

 **Example Command**
.. code-block::

    prokka --outdir {params.outdir} --locustag {params.locustag} \
    --compliant --prefix {params.locustag} {input.scaffolds} --force


-----------------------
Alternative Approach
-----------------------
Alternatively


--------------------------
Assessing Assembly Quality
--------------------------



