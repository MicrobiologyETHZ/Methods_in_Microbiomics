
# Metagenome Assembly

Technical advances in sequencing technologies in recent decades have allowed detailed investigation of complex microbial communities without the need for cultivation. Sequencing of microbial DNA extracted directly from environmental or host-associated samples have provided key information on microbial community composition. These studies have also allowed gene-level characterization of microbiomes as the first step to understanding the communities' functional potential. Furthermore, algorithmic improvements, as well as increased availability of computational resources, make it now possible to reconstruct whole genomes from metagenomic samples (metagenome-assembled genomes (MAGs)). Methods for microbial community composition analysis are discussed in [here](../profiling/metagenomes). Here we describe building [metagenomic assemblies](massembly), as well as building [gene catalogs](gene_catalogs), and [MAGs](mags) from metagenomic data.

(massembly)= 
## Metagenomic Assembly

```{mermaid}

    flowchart LR
            id1( Metagenomic<br/>Assembly) --> id2(data preprocessing)
            id2 --> id3(metagenomic assembly<br/>fa:fa-cog metaSPAdes)
            id3 --> id4(gene catalogs)
            id3 --> id5(MAGs)
            classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
            style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
            style id2 fill:#F78A4A,stroke:#F8F7F7,stroke-width:1px
            class id3,id4,id5 tool

```

1. **Data Preprocessing**. Before proceeding to the assembly, it is important to preprocess the raw sequencing data. Standard preprocessing protocols are described in [Data Preprocessing](../preprocessing/preprocessing). In addition to standard quality control and adapter trimming, we also suggest normalization with **bbnorm.sh** and merging of paired-end reads (see [Data Preprocessing](../preprocessing/preprocessing) for more details).

2. **Metagenomic Assembly**. Following data preprocessing, we use clean reads to perform a metagenomic assembly using **metaSPAdes**. metaSPAdes is part of the [SPAdes](https://github.com/ablab/spades) assembly toolkit. 
**Example Command**:

```bash
    mkdir metag_assembly
    for i in 1 2 3
      do
        mkdir metag_assembly/metag$i
        metaspades.py -t 4 -m 10 --only-assembler \
        --pe1-1 reads/metag$i.1.fq.gz \
        --pe1-2 reads/metag$i.2.fq.gz \
        -o metag_assembly/metag$i
      done
```

|                      |                                                                                                          |
|:---------------------|:---------------------------------------------------------------------------------------------------------|
| ``-t``               | Number of threads                                                                                        |
| ``-m``               | Set memory limit in Gb; spades will terminate if that limit is reached                                   |
| ``--only-assembler`` | Run assembly module only (spades can also perform read error correction, <br> this step will be skipped) |
| ``--pe1-1``          | Forward reads                                                                                            |
| ``--pe1-2``          | Reverse reads                                                                                            |
| ``--pe1-s``          | Unpaired reads                                                                                           |
| ``--pe1-m``          | Merged reads                                                                                             |
| ``-o``               | Specify output directory                                                                                 |

```{note}

**Computational Resources** needed for metagenomic assembly will vary significantly between datasets. In general, metagenomic assembly requires a lot of memory (usually > 100 Gb). You can use multiple threads (16-32) to speed up the assembly. Because test data set provided is very small, merging of the pair-end reads was not necessary (see :doc:`../preprocessing/preprocessing`). It is helpful when working with real data - don't forget to include the merged and singlton files with ``--pe1-m`` and ``--pe1-s`` options.
```


3. **Filtering**. Following the assembly, we generate some assembly statistics using **assembly-stats**, and filter out contigs that are < 1 kbp in length. The script we use for scaffold filtering can be found here: {download}`scaffold_filter.py <../scripts/scaffold_filter.py>`. It is also included in the test dataset for this section.

**Example Command: filtering**

```bash
    cd metag_assembly
      for i in 1 2 3
        do
          python ../scaffold_filter.py metag$i scaffolds metag$i/scaffolds.fasta metag$i META
        done

```

|                            |                                                                          |
|:---------------------------|:-------------------------------------------------------------------------|
| ``metag1``                 | Sample name                                                              |
| ``scaffolds``              | Sequence type (can be contigs, <br> scaffolds or transcripts)            |
| ``metag1/scaffolds.fasta`` | Input assembly to filter                                                 |
| ``metag1``                 | Prefix for the output file                                               |
| ``META``                   | Type of assembly (META for metagenomics <br> or ISO for isolate genomes) |


**Example Command: stats**

```bash
      for i in 1 2 3
        do
          assembly-stats -l 500 -t <(cat metag$i/metag$i.scaffolds.min500.fasta) \
          > metag$i/metag$i.assembly.stats
        done


```
|          |                                         |
|:---------|:----------------------------------------|
| ``-l``   | Minimum length cutoff for each sequence |
| ``-t``   |        Print tab-delimited output       |


4. The metagenomic scaffolds generated in step 2 can now be used to build and/or profile [gene catalogs](gene_catalogs) or to construct [MAGs](mags).


(gene_catalogs)=
## Gene Catalogs

Gene catalog generation and profiling (i.e. gene abundance estimation) can provide important insights into the community's structure, diversity and functional potential. This analysis could also identify relationships between genetic composition and environmental factors, as well as disease associations.

```{note}

Integrated catalogs of reference genes have been generated for many ecosystems (e.g. [ocean](https://doi.org/10.1016/j.cell.2019.10.014), [human gut](https://doi.org/10.1038/s41587-020-0603-3), and [many others](https://doi.org/10.1038/s41586-021-04233-4)) and might be a good starting point for the analysis.

```

### Building

This protocol will allow you to create a de novo gene catalog from your metagenomic samples.

```{mermaid}
    flowchart LR
        id1( Building a<br/>Gene Catalog) ---> id2(gene calling<br/>fa:fa-cog prodigal)
        id2 ---> id3(gene dereplication<br/>fa:fa-cog CD-HIT)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        class id2,id3 tool

```

1. **Gene calling**. We use [**prodigal**](https://github.com/hyattpd/Prodigal) to extract protein-coding genes from metagenomic assemblies (using **scaffolds** >= 500 bp as input). Prodigal has different gene prediction modes with single genome mode as default. To run prodigal on metagenomic data, we add the ``-p meta`` option. This will produce a fasta file with amino acid sequences (.faa), nucleotide sequences (.fna) for each gene, as well as an annotation file (.gff).

**Example Command**

```bash
  for i in 1 2 3
          do
            prodigal -a metag$i/metag$i.faa -d metag$i/metag$i.fna -f gff \
            -o metag$i/metag$i.gff -c -q -p meta \
            -i metag$i/metag$i.scaffolds.min500.fasta
          done


```

``-a``           Specify protein translations file
``-d``           Specify nucleotide sequences file
``-f``           Specify output format: gbk: Genbank-like format (Default); gff: GFF format; sqn: Sequin feature table format; sco: Simple coordinate output
``-o``           Specify output file, default stdout
``-c``           Closed ends, do not allow partial genes at edges of sequence
``-q``           Run quietly (suppress logging output)
``-p``           Specify mode: single or meta
``-i``           Input FASTA or Genbank file


(mags)= 
## MAGs

```{admonition} Here's my title
:class: warning

Here's my admonition content
```


```{mermaid}

   flowchart LR
        id1( Metagenomic<br/>Assembly) --> id2(data preprocessing)
        id2 --> id3(metagenomic assembly<br/>fa:fa-cog metaSPAdes)
        id3 --> id4(gene catalogs)
        id3 --> id5(MAGs)
        classDef tool fill:#96D2E7,stroke:#F8F7F7,stroke-width:1px;
        style id1 fill:#5A729A,stroke:#F8F7F7,stroke-width:1px,color:#fff
        style id2 fill:#F78A4A,stroke:#F8F7F7,stroke-width:1px
        class id3,id4,id5 tool

```

```bash
mkdir metag_assembly
for i in 1 2 3
  do
    mkdir metag_assembly/metag$i
    metaspades.py -t 4 -m 10 --only-assembler \
    --pe1-1 reads/metag$i.1.fq.gz \
    --pe1-2 reads/metag$i.2.fq.gz \
    -o metag_assembly/metag$i
  done
```

