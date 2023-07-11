======================
Gene catalogue profiling (Metatranscriptomics & Metagenomics)
======================

Protocol provided by Guillem Salazar.

-----------------------------------------
Dataset
-----------------------------------------
In this documentation we will combine metatranscriptomics with metagenomics by using a genome as a reference for the transcriptome.

To analyze the data, we first need to remove the following sources of bias:

* Differences in gene length between genes.
* Differences in sequencing depth between samples.
* Differences in genome size distribution between samples.
* Compositionality: The number of inserts for a given gene in a given sample is in itself arbitrary and can only be interpreted relative to the rest of the genes in the sample.

We normalize using the following steps:

1. Divide the insert counts by the gene length for each gene in each sample.
2. Compute the total insert count of 10 universal single-copy marker genes (MGs) in each sample.
3. Compute the median across the 10 MGs in each sample.
4. Divide the gene-length normalized abundances by this median in each sample.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Why do we divide by the 10 Marker Genes?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The 10 Marker Genes are:

* Universal: present in "all" prokaryotes
* Single-copy: always present once per cell (genome)
* Long enough for a satisfactory profiling
* Are housekeeping genes

The marker genes serve as a gene catalogue for our analysis. A gene catalogue is a data structure that organizes genes and provides a reference for standardized analysis of the microbiomes across samples and studies. Further information can be found in :doc:`../assembly/metagenomic_workflows`.

Therefore, their abundance should be well correlated and the median of their abundances should correlate with the sequencing depth.

The median of their abundances is a good proxy for the amount of cells captured in a given metagenome. The per-cell normalization should account for differences in genome sizes between samples and the per-cell normalization is in fact a way for controlling for composition. The result has a biologically meaningful unit: *gene copies per total cell in the community*.

The dataset used in this tutorial is from the article `Gene Expression Changes and Community Turnover Differentially Shape the Global Ocean Metatranscriptome, Salazar et al <https://doi.org/10.1016/j.cell.2019.10.014>`_.
The data can be downloaded :download:`here <../downloads/metatranscriptomics_tutorial.zip>`.


-----------------------------------------
Performing normalization based on gene length and the abundance of the 10 MGs
-----------------------------------------

This script normalizes the data based on gene length and abundance of 10 marker genes. In this example we use the following marker genes:

.. image:: ../images/metatranscriptomics_marker_genes.jpg

.. note::
    Running this code with the given data needs a significant amount of resources. We recommend running it on a server.

.. code-block:: r

    library(data.table)
    library(tidyverse)
    library(R.utils)

    # Define the KOs corresponding to the 10 MGs
    mgs<-c("K06942", "K01889", "K01887", "K01875", "K01883", "K01869", "K01873", "K01409", "K03106", "K03110")

    # Load the gene catalog (insert counts)
    gc_profile<-fread("gene_profile_tara.tsv.gz",sep="\t",header=T,data.table = F,tmpdir=".")
    sample_info<-fread("sample_info_tara.tsv",sep="\t",header=T,data.table = F,tmpdir=".")

    # Assign the median gene length to the -1 (the category designing all reads not mapping to any gene)
    if (length(which(gc_profile$length<0))>0){
      gc_profile$length[which(gc_profile$length<0)]<-median(gc_profile$length[which(gc_profile$length>0)])
    }

    # Build a gene-length normalized profile
    gc_profile_lengthnorm<-gc_profile[,1:4]
    for (i in 5:ncol(gc_profile)){
      cat("Normalizing by gene length: sample",colnames(gc_profile)[i],"\n")
      tmp<-gc_profile[,i]/gc_profile$length %>%
        as.data.frame()
      colnames(tmp)<-colnames(gc_profile)[i]
      gc_profile_lengthnorm<-gc_profile_lengthnorm %>%
        bind_cols(tmp)
    }

    # Build a MGs normalized profile
    gc_profile_lengthnorm_mgnorm<-gc_profile_lengthnorm[,1:4]
    for (i in 5:ncol(gc_profile_lengthnorm)){
      cat("Normalizing by 10 MGs: sample",colnames(gc_profile_lengthnorm)[i],"\n")
      mg_median<-gc_profile_lengthnorm %>%
        select(KO,abundance=all_of(colnames(gc_profile_lengthnorm)[i])) %>%
        filter(KO %in% mgs) %>%
        group_by(KO) %>% summarise(abundance=sum(abundance)) %>%
        ungroup() %>% summarise(mg_median=median(abundance)) %>%
        pull()
      tmp<-gc_profile_lengthnorm[,i]/mg_median
      tmp<-tmp %>% as.data.frame()
      colnames(tmp)<-colnames(gc_profile_lengthnorm)[i]
      gc_profile_lengthnorm_mgnorm<-gc_profile_lengthnorm_mgnorm %>%
        bind_cols(tmp)
    }

    # Save profiles and compress
    fwrite(gc_profile_lengthnorm,"gene_profile_tara_lengthnorm.tsv",sep="\t")
    gzip("gene_profile_tara_lengthnorm.tsv")
    fwrite(gc_profile_lengthnorm_mgnorm,"gene_profile_tara_lengthnorm_percell.tsv",sep="\t")
    gzip("gene_profile_tara_lengthnorm_percell.tsv")

-----------------------------------------
Showing the effect of the normalization
-----------------------------------------
Here, we visualize the effect of the normalization based on length and abundance of marker genes. Using this script we create the following plots:

.. image:: ../images/K03040_K03043_comparison.jpg
.. image:: ../images/mgs_vs_seqdepth.jpg
.. image:: ../images/mgs_pairwise_corr.jpg

.. code-block:: r

    library(data.table)
    library(tidyverse)
    library(patchwork)
    library(GGally)
    library(R.utils)

    # Define the KOs corresponding to the 10 MGs
    mgs<-c("K06942", "K01889", "K01887", "K01875", "K01883", "K01869", "K01873", "K01409", "K03106", "K03110")

    # Load the raw count and gene length normalized profiles and the sample information
    gc_profile<-fread("gene_profile_tara.tsv.gz",sep="\t",header=T,data.table = F,tmpdir=".")
    gc_profile_lengthnorm<-fread("gene_profile_tara_lengthnorm.tsv.gz",sep="\t",header=T,data.table = F,tmpdir=".")
    sample_info<-fread("sample_info_tara.tsv",sep="\t",header=T,data.table = F,tmpdir=".")

    # Compute the abundance of K03040 and K03043 with and without gene-length normalization
    rp_ab<-gc_profile %>%
      select(-reference,-length,-Description) %>%
      filter(KO %in% c("K03040","K03043")) %>%
      pivot_longer(-KO,names_to = "sample",values_to = "inserts") %>%
      filter(sample %in% sample_info$sample_metag) %>%
      group_by(KO,sample) %>% summarise(inserts=sum(inserts)) %>%
      pivot_wider(names_from = "KO",values_from = "inserts")

    rp_ab_lengthnorm<-gc_profile_lengthnorm %>%
      select(-reference,-length,-Description) %>%
      filter(KO %in% c("K03040","K03043")) %>%
      pivot_longer(-KO,names_to = "sample",values_to = "inserts_lengthnorm") %>%
      filter(sample %in% sample_info$sample_metag) %>%
      group_by(KO,sample) %>% summarise(inserts_lengthnorm=sum(inserts_lengthnorm)) %>%
      pivot_wider(names_from = "KO",values_from = "inserts_lengthnorm")

    g1<-ggplot(data=rp_ab,aes(x=K03040,y=K03043)) +
      geom_point(alpha=0.5) +
      geom_abline(slope = (1342/329)) +
      geom_abline(linetype=2) +
      xlim(range(rp_ab$K03040,rp_ab$K03043)) +
      ylim(range(rp_ab$K03040,rp_ab$K03043)) +
      xlab("K03040: rpoA\n(DNA-directed RNA polymerase subunit alpha)") +
      ylab("K03043: rpoB\n(DNA-directed RNA polymerase subunit beta)") +
      labs(title="Insert counts",subtitle="Slope ~ 4 which corresponds to the ratio of gene lengths\n(K03040: 1,342 aa; K03043: 329 aa in E. coli K-12)") +
      coord_fixed() +
      theme_bw() +
      theme(plot.subtitle = element_text(size=7))
    g2<-ggplot(data=rp_ab_lengthnorm,aes(x=K03040,y=K03043)) +
      geom_point(alpha=0.5) +
      geom_abline(linetype=2) +
      xlim(range(rp_ab_lengthnorm$K03040,rp_ab_lengthnorm$K03043)) +
      ylim(range(rp_ab_lengthnorm$K03040,rp_ab_lengthnorm$K03043)) +
      xlab("K03040: rpoA\n(DNA-directed RNA polymerase subunit alpha)") +
      ylab("K03043: rpoB\n(DNA-directed RNA polymerase subunit beta)") +
      labs(title="Gene-length normalized insert counts",subtitle="Slope ~ 1 once insert counts are corrected for differences\nin gene lengths") +
      coord_fixed() +
      theme_bw() +
      theme(plot.subtitle = element_text(size=7))
    g<-g1 | g2
    ggsave("K03040_K03043_comparison.pdf",g,width=unit(10,"cm"),height=unit(4.5,"cm"))

    # Compute the abundance of the 10MGs and correlate to sequencing depth
    mgs_ab_lengthnorm<-gc_profile_lengthnorm %>%
      select(-reference,-Description,-length) %>%
      filter(KO %in% mgs) %>%
      pivot_longer(-KO,names_to = "sample",values_to = "inserts_lengthnorm") %>%
      group_by(KO,sample) %>% summarise(inserts_lengthnorm=sum(inserts_lengthnorm)) %>%
      ungroup() %>% group_by(sample) %>% summarise(median_mgs=median(inserts_lengthnorm)) %>%
      inner_join(sample_info,by=c("sample"="sample_metag"))


    g3<-ggplot(data=mgs_ab_lengthnorm,aes(x=sample_metag_nreads,y=median_mgs)) +
      geom_point(alpha=0.5) +
      #geom_smooth(method = "lm") +
      #scale_x_log10() +
      #scale_y_log10() +
      xlab("Sequencing depth (number of reads)") +
      ylab("Median abundance of the 10 universal\nand single-copy marker genes") +
      theme_bw() +
      theme(legend.title = element_blank())

    ggsave("mgs_vs_seqdepth.pdf",g3,width=unit(7,"cm"),height=unit(4,"cm"))

    # Compute the abundance of the 10MGs and their autocorrelation
    mgs_ab_lengthnorm<-gc_profile_lengthnorm %>%
      select(-reference,-Description,-length) %>%
      filter(KO %in% mgs) %>%
      pivot_longer(-KO,names_to = "sample",values_to = "inserts_lengthnorm") %>%
      group_by(KO,sample) %>% summarise(inserts_lengthnorm=sum(inserts_lengthnorm)) %>%
      inner_join(sample_info,by=c("sample"="sample_metag")) %>%
      select(KO,sample,inserts_lengthnorm) %>%
      pivot_wider(names_from = "KO",values_from = "inserts_lengthnorm")

    g4<-ggpairs(data=mgs_ab_lengthnorm %>% column_to_rownames("sample")) +
      scale_x_log10() +
      scale_y_log10()

    ggsave("mgs_pairwise_corr.pdf",g4,width=unit(10,"cm"),height=unit(10,"cm"))



-----------------------------------------
Combining Metatranscriptomic and Metagenomic Data
-----------------------------------------
In this section we combine metatranscriptomic and metagenomic data and create the following plot:

.. image:: ../images/K03704.jpg

.. code-block:: r

    library(data.table)
    library(tidyverse)
    library(patchwork)
    library(GGally)
    library(R.utils)

    # Load normalized profile
    gc_profile<-fread("gene_profile_tara_lengthnorm_percell.tsv.gz",sep="\t",header=T,data.table = F,tmpdir=".")
    sample_info<-fread("sample_info_tara.tsv",sep="\t",header=T,data.table = F,tmpdir=".")

    # Build a KO profile by adding up all genes with the same KO annotation
    ko_profile<-gc_profile %>%
      group_by(KO) %>% summarise(across(starts_with("TARA"),sum)) %>%
      as.data.frame()

    # Compute the gene abundance, transcript abundance and expression for the pairs of metaG-metaT samples
    # The expression is just the ratio of transcript_abundance to gene_abundance
    tmp_sample_info<-sample_info %>%
      select(sample_metag,sample_metat) %>%
      mutate(sample_pair=paste(sample_metag,sample_metat,sep="-"))
    tmp_metag<-ko_profile %>%
      select(KO,all_of(tmp_sample_info$sample_metag)) %>%
      pivot_longer(-KO,names_to = "sample_metag",values_to = "gene_abundance")
    tmp_metat<-ko_profile %>%
      select(KO,all_of(tmp_sample_info$sample_metat)) %>%
      pivot_longer(-KO,names_to = "sample_metat",values_to = "transcript_abundance")
    final_profile<-tmp_sample_info %>%
      left_join(tmp_metag,by="sample_metag") %>%
      left_join(tmp_metat,by=c("KO","sample_metat")) %>%
      mutate(expression=transcript_abundance/gene_abundance)

    # Plot the gene abundance, expression and transcript abundance of K03704: cspA: cold shock protein
    toplot<-final_profile %>%
      filter(KO=="K03704") %>%
      left_join(sample_info,by=c("sample_metag","sample_metat"))

    g_metat<-ggplot(data=toplot,aes(y=transcript_abundance,x=`temperature [°C]`,color=depth_layer)) +
      geom_point() +
      geom_smooth(method = "gam",se = T) +
      #scale_y_log10() +
      #coord_flip() +
      scale_color_manual(values=c("darkgreen","darkblue")) +
      ylab("Transcript abundance") +
      theme_bw() +
      theme(legend.position = "none")
    g_metag<-ggplot(data=toplot,aes(y=gene_abundance,x=`temperature [°C]`,color=depth_layer)) +
      geom_point() +
      geom_smooth(method = "gam",se = T) +
      #scale_y_log10() +
      #coord_flip() +
      scale_color_manual(values=c("darkgreen","darkblue")) +
      ylab("Gene abundance") +
      theme_bw() +
      theme(legend.position = "none")
    g_exp<-ggplot(data=toplot,aes(y=expression,x=`temperature [°C]`,color=depth_layer)) +
      geom_point() +
      geom_smooth(method = "gam",se = T) +
      #scale_y_log10() +
      #coord_flip() +
      scale_color_manual(values=c("darkgreen","darkblue")) +
      ylab("Gene expression") +
      theme_bw() +
      theme(legend.position = "top",legend.title = element_blank())
    g<-g_metag | g_exp | g_metat
    ggsave("K03704.pdf",g)
