==================================
Taxonomic Profiling of Metagenomes
==================================

Taxonomic profiling of complex microbial communities is an essential first step in the  investigation of relationship between community composition and environmental and/or health factors. The most common approach to community profiling is amplification and classificaton the 16S rRNA gene. Methods related to 16S rRNA analysis are discussed in detail in :doc:`../profiling/16S`. Recently shotgun metagenomic sequencing has started to replace the amplicon based approaches, as it provides higher resolution information about the microbial community, and resolves some of the biases associated with 16S approach. A number of software tools have been developed to taxonomically profile metagenomic samples. These tools have been profiled in detail in recent studies<add link>. Here we're going to talk about the use of :ref:`mOTUs` and :ref:`mTAGs` for taxonomic profiling.

--------
mOTUs
--------

- Overview
- Database

.. note::

    Sample data and conda environment file for this section can be found :download:`here <../downloads/metag.tar.gz>`. See the :ref:`tutorials` section for instructions on how to unpack the data and create the conda environment. ...


.. mermaid::

   flowchart LR
        id1( Taxonomic profiling<br/>with mOTUs) --> id2(data preprocessing)
        id2 --> id3(profile<br/>fa:fa-cog mOTUs)
        id3 --> id4(merge<br/>fa:fa-cog mOTUs)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        style id2 fill:#F78A4A,stroke:#F8F7F7,stroke-width:1px
        class id3,id4 tool


1. **Data Preprocessing**. Before proceeding to the assembly, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in :doc:`../preprocessing/preprocessing`. In addition to standard quality control and adapter trimming, we also suggest merging of paired-end reads (see :doc:`../preprocessing/preprocessing` for more details). Using merged reads increases speed and accuracy.

2. **Profile**.

.. code-block:: bash

    motus profile -f {in.1.fq.gz } -r {in.2.fq.gz} -s {in.s.fq.gz},{in.m.fq.gz} -n {sample} \
-o {out.motus} -I {out.bam} -t {threads} -v 100 -y {counting_method} -c


==============    =====================================================================================================
``-f``
``-r``
``-s``
``-n``
``-o``
``-I``
``-t``
``-v``
``-y``
``-c``
``-k``
==============    =====================================================================================================

3. **Merge**.

.. code-block:: bash

    motus merge -i ERR479298s-default.motus,bp1_precomputed/ERR479045- default.motus -o merged.motus


--------
mTAGs
--------

