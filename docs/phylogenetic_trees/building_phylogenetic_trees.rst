============================
Building phylogenetic Trees
============================

Protocol provided by Marija Dmitrijeva.

Phylogenetic trees are a means to explore evolutionary relationships. For example, they can be used to gain insight
into major events such as `the emergence of eukaryotic cells`_,
establish `taxonomic relationships for classification of novel taxa`_,
or detect `signatures of horizontal gene transfer`_. In the phylogenetic tree, leaves correspond to
present-day biological entities, and internal nodes correspond to their ancestors prior to speciation, i.e.,
prior to one lineage splitting into two divergent lineages.

.. _the emergence of eukaryotic cells: https://www.doi.org/10.1038/s41586-023-06186-2

.. _taxonomic relationships for classification of novel taxa: https://www.doi.org/10.1038/nbt.4229

.. _signatures of horizontal gene transfer: https://www.doi.org/10.21203/rs.3.rs-3062985/v1

A key concept in phylogenetic tree construction is sequence homology, i.e., the sequences used to build the
phylogenetic tree must share common ancestry. Divergent sequences sharing ancestry could stem from various
evolutionary events: speciation (resulting in orthologs), gene duplication (resulting in paralogs), or
horizontal gene transfer (resulting in xenologs). Often, high nucleotide or amino acid sequence similarity
can be used to infer homology.

The steps to build a phylogenetic tree are the following:

.. mermaid::

   flowchart LR
        id1(Homologous<br/>Sequence<br/>Identification) --> id2("Multiple<br/>Sequence<br/>Alignment<br/>(MAFFT or MUSCLE)")
        id2 --> id3("Phylogenetic<br/>Inference<br/>(FastTree or RaxML NG)")
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        class id1,id2,id3 tool

-----------------------------------
Homologous Sequence Identification
-----------------------------------





