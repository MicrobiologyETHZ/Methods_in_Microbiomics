=================================
SNV Analysis on Metagenomic Data
=================================

Protocol provided by Aiswarya Prasad.

Single genomes and metagenome-assembled genomes (MAGs) provide a single snapshot in time and space of a community of
bacteria. They can also be thought of as an average representation of a group of very closely related microbes (strains
of a species) found in that sample or pool of samples. Profiling single-nucleotide variants (SNVs) provides deeper
insight into the sub-populations of the community. With sequencing becoming cheaper and faster, it is now possible
to sequence metagenomic samples deep enough to recover MAGs and use those to detect SNVs in those samples.

SNV profiling is cheaper and more straightforward for simple microbial communities with few species or genera. A crude
way of assessing the minimum depth (e.g., 10x) needed a priori is to consider the number of species/genomes you are
interested in. Suppose you are interested in **10** genomes where the length of these genomes is **2** Mb with reads of
**150** bp, you would calculate the following:

.. math::

    2 \: Mb * 10 \: genomes * 10 \: coverage \: / \: 150 \: read \: length \sim 1.3 \: million \: reads \: (\sim 200 \: Mb \: of \: data)

In this case, the above calculation underestimates the coverage, because it assumes that the ten genomes
have similar abundance, which is rarely true. As a result, more abundant genomes will have more than 10x coverage, while
the low-abundance genomes of interest will have little to none. Therefore, it is best to sequence deeper
than this calculation suggests. In addition, if you have pooled samples, you must sequence deeper than you would
sequence individual samples, and you would also lose information about strain-level variation between individuals.

.. note:: 
    **What are bacterial species?**

    While the boundaries specifying units of diversity in microbial communities are not easily defined, bacterial genomes can still be categorized into units based on their genome sequence. One popular
    method to do this is to define clusters of genomes that share >95% `average nucleotide identity (ANI)`_. These clusters in most
    cases correspond to genomes of the same species. This approach is becoming more prevalent, especially due to cheaper
    sequencing methods that allow us to recover MAGs from complex samples, as well as tools that allow for pairwise comparisons
    of genomes.

.. _average nucleotide identity (ANI): https://doi.org/10.1038/s41467-018-07641-9

In this protocol, we describe an approach to analyze population-level diversity using `inStrain`_, which is documented
in detail `here`_.

.. _inStrain: https://doi.org/10.1038/s41587-020-00797-0

.. _here: https://instrain.readthedocs.io/en/latest/index.html

----------------------
Aim
----------------------
The goal of this analysis is to profile SNVs in a metagenomic sample in a database-independent manner, starting from bam
files obtained by mapping trimmed reads from each sample to a database of MAGs recovered from all the samples
de-replicated into 95% clusters (recommended). For guidelines on how to preprocess sequencing data and assemble MAGs, please refer to :doc:`../preprocessing/preprocessing`
and :doc:`../assembly/metagenomic_workflows`.


----------------------
Overview
----------------------
We must complete a few steps to prepare the data for `inStrain`. Including recovering the set of MAGs from the samples of
interest, dereplicating the MAGs, and gene calling or annotation of the MAGs (e.g., using Prokka). In addition, it is
helpful to create a metadata file containing information about the MAGs, such as which dereplicated group they belong
to, their quality score, etc. The tutorial :doc:`../assembly/metagenomic_workflows` covers building metagenomic assembly.

`inStrain` outlines the considerations to be made when preparing the input in `their tutorial on establishing and
evaluating genome databases`_. This is a very useful place to gain an understanding of some important concepts and
the reasons behind choosing this approach. This tutorial includes examples of code that can be used to achieve this.

The tutorials on the `inStrain documentation`_ page are also very helpful. This tutorial presents an alternate scenario
parallel to Tutorial #2 on that page.

.. _their tutorial on establishing and evaluating genome databases: https://instrain.readthedocs.io/en/latest/important_concepts.html?highlight=drep#establishing-and-evaluating-genome-databases

.. _inStrain documentation: https://instrain.readthedocs.io/en/latest/tutorial.html#tutorials

Create a Database of Metagenome-assembled Genomes (MAGs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this step, we create a database of MAGs containing one MAG to represent each species to infer SNVs. This can be
done by choosing the highest quality MAG from each cluster of MAGs that are on average, e.g. 95% similar to each other.
The MAGs can be clustered using the tool `dRep`_. Once you have collected and renamed all medium and good quality (>50%
completeness and <5% contamination - as a rule of thumb) MAGs, make a tab-separated file containing the
columns named:

- :code:`Bin Id` - genome name with the extension
- :code:`Completeness` - from `checkM`
- :code:`Contamination` - from `checkM`

You can do this using a script that parses the checkM output for the list of MAGs of at least medium
quality.

.. _dRep: https://drep.readthedocs.io/en/latest/index.html

.. code-block:: bash

    dRep dereplicate /your/output/directory -g ./folder/file.fa -comp 0 -con 1000 --clusterAlg average \
            --genomeInfo drep_genome_info.tsv -sa 0.95 -nc 0.2 -p 1 --debug

The result that we are interested in will be found in :code:`/your/output/directory/data_tables/`

- :code:`Sdb.csv` lists each MAG and its score.
- :code:`Cdb.csv` contains the information about which cluster each MAG ends up in.

Using this information, collect the best scoring MAG for each dRep cluster (on the second column of :code:`Cdb.csv` named
:code:`secondary_cluster`) and concatenate them into the file :code:`MAG_rep_database.fa`. Ensure that the scaffolds have unique
names across MAGs. In the rep database fasta, you will simply put all the scaffolds with each having its own header
and not necessarily the information about which MAG it is from. This information will be separately summarized in
another file for later steps.

Annotate the representative Database (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you would like to have gene-level information in inStrain you can include a file specifying the positions of genes
in the MAGs obtained using prodigal (this is optional for `inStrain`).

.. code-block:: bash

    prodigal -i MAG_rep_database.fa -d MAG_rep_database.fna -a MAG_rep_database.faa -p meta &> my_log_file.log

Make Scaffolds to bin File
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: python

    with open(output.scaffold_to_bin_file, "w") as f:
        for mag in input.rep_mags:
            with open(mag, "r") as m:
                for line in m:
                    mag_name = os.path.basename(mag).split(".")[0]
                    if line.startswith(">"):
                        scaffold = line.strip().split(">")[1]
                        f.write(f"{scaffold}\t{mag_name}\n")

Map reads to MAG Database
^^^^^^^^^^^^^^^^^^^^^^^^^
Ensure that you use `bowtie2` for this, as recommended by `inStrain`. Avoid BWA (even though it might be your favorite
aligner) as `inStrain` may have issues handling the way that it calculates insert size, and the BWA documentation is
unclear about how this is performed.

.. code-block:: bash

    bowtie2-build mag_rep_database.fa mag_rep_database.fa &> bowtie2_build.log # bowtie index
    # map to rep MAGs
    bowtie2 -X 1000 -x mag_rep_database.fa -1 sample_R1_repaired.fastq.gz -2 sample_R2_repaired.fastq.gz | samtools view -bh - | samtools sort - > sample_bowtie.bam
    samtools flagstat sample_bowtie.bam > sample_bowtie_flagstat.tsv

Make the inStrain Profile
^^^^^^^^^^^^^^^^^^^^^^^^^
Output and parameter information is well-documented in inStrain - run using db_mode if you wish to run `inStrain compare`
later. This makes it much faster.

.. code-block:: bash

    inStrain profile sample_bowtie.bam mag_rep_database.bam -o /your/output/directory/ -p 8 -g mag_rep_database_genes.fna \
        --max_insert_relative 5 -s scaffold_to_bin_file.tsv
    inStrain plot -i inStrain_profile_object -pl a -p 16
    inStrain profile sample_bowtie.bam mag_rep_database.fa -o /path/to/output/folder/sample -p 8 \
        -g mag_rep_database_genes.fna --max_insert_relative 5 --database_mode -s scaffold_to_bin_file.tsv

Run `inStrain compare`
^^^^^^^^^^^^^^^^^^^^^^^^^
This can be run on all profiles together, especially if you did not have a lot of samples, but for datasets including a
large number of samples, it will be more efficient to run this in parallel for each species at a time by using the
:code:`--genome` option to specify one genome at a time.

.. code-block:: bash

    inStrain compare -i inStrain_profile -s scaffold_to_bin_file.tsv -p 8 -o /your/output/directory/ \
                --database_mode --genome mag.stb



