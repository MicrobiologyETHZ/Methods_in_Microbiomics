==================================================
Old guide
==================================================

-----------------------------------------
The datasets
-----------------------------------------
RNA was extracted, processed and analyzed with a hybrid approach where both de-novo assembly of the data is used together with the use of the genome as reference.
Out of this procedure, transcript abundance tables are produced where the number of transcripts is quantified for each gene in the sample. As the main goal of the present guide is to analyze the response of the native soil community, all genes from the inoculum were further excluded.
Genes are finally assigned a metabolic function by using the KEGG database. All genes belonging to the same function are added up and thus the gene abundance table is then transformed into a functional abundance table. For example, the transcript abundance of all genes from different organisms that encode the pyruvate kinase (designated as K00873 in the KEGG database) are added up. In this way we quantify the abundance of each metabolic function in each sample (regardless of the microbe encoding for it).
This functional abundance table and an additional table with the sample information are the two main ingredients for the following analysis.


-----------------------------------------
Loading and formatting the data
-----------------------------------------
The next thing we need to do is to load the data, specifically:

#.   The functional abundance table. We will load it and save it as the object kotab.
#.   The table with the information on the samples. We will load it and save it as the object metadata.
#.   And an additional table with the description of each function.We will load it and save it as the object kegg_annot.

.. code-block:: r

    kotab<-fread("kotab_insertcounts.tsv",sep="\t",header=T,data.table = F)

    metadata<-fread("metadata.tsv",sep="\t",header=T,data.table = F)

    kegg_annot<-fread("kegg-annotations.tsv",sep="\t",header=T,data.table = F) %>%
        distinct(KO,Description)

Once the data is loaded, we need to proceed with some basic transformation of the data. Specifically, the functional abundance table needs to be converted into a numeric matrix after moving the sample identifier to the rownames. For that, we will create a new object (kotab_matrix):

.. code-block:: r

    kotab_matrix<-kotab %>%
        column_to_rownames(var="replicate_name") %>%
        as.matrix()

-----------------------------------------
Ordination
-----------------------------------------
The primary aim of ordination is to represent multiple objects in a reduced number of orthogonal (i.e., independent) axes. Ordination plots are particularly useful for visualizing the similarity among objects. For example, in the context of the current dataset, samples that are closer in ordination space have functional compositions that are more similar to one another than samples that are further apart in ordination space.

There are various ordination techniques that can be applied to multivariate microbiome data. Common methods include: Principal Components Analysis (PCA), Correspondence Analysis (CA), Principal Coordinates Analysis (PCoA), Factor Analysis (FA), and Nonmetric Multidimensional Scaling (NMDS). The distinction between these is not the goal of this guide.

Principal components analysis (PCA) is one of the oldest ordination techniques, and probably one that you may already know. It provides graphs that show the Euclidean distance between sites. No other distances can be investigated with PCA. This ordination method is not ideal for analysis of information on species, gene or transcript abundances because of the limitations of the Euclidean distance with this type of data.

While only euclidean distances may be used with PCA, any kind of dissimilarity measure is suitable for the Principal Coordinates Analysis (PCoA). Its use is similar to the use of PCA, however, it can handle non-euclidean distances as the one we need here.

The first thing we need to do is to compute how similar each pair of sample is based on their functional compositions. We will compute the Bray-Curtis dissimilarity index. This will create a matrix (that we will store as kotab_bc) with a value between 0 and 1 for each pair of samples. The closer the value is to 0 the more similar are two samples, the closer it is to 1 the more different they are.

.. code-block:: r

    kotab_bc<-vegdist(kotab_matrix)

This dissimilarity matrix can be visualized through the PCoA. We will run the PCoA analysis and plot all samples along the first two axes.

.. code-block:: r

    pcoa_kotab<-pcoa(kotab_bc)

    toplot_kotab<-data.frame(pcoa_kotab$vectors) %>%
      select(Axis.1:Axis.5) %>%
      rownames_to_column(var="replicate_name") %>%
      left_join(metadata,by="replicate_name")

    ggplot(data=toplot_kotab,aes(x=Axis.1,y=Axis.2,shape=variable)) +
      geom_point(alpha=0.7,size=4) +
      theme_bw() +
      xlab(paste("Axis 1 ","(",round(100*pcoa_kotab$values$Relative_eig[1],1),"%)",sep="")) +
      ylab(paste("Axis 2 ","(",round(100*pcoa_kotab$values$Relative_eig[2],1),"%)",sep=""))


-----------------------------------------
Hypothesis testing
-----------------------------------------
While the ordination techniques above mentioned are very useful and are the logical first step to explore high-dimensional data, no statistically supported conclusions can be obtained from them. Ordination techniques often are supplemented by hypothesis testing techniques (i.e. statistical tests). The main hypothesis to be tested with this kind of data is wether the similarity between samples is organized in pre-defined groups. That is, if different clusters of samples exists based on their similarity.

Permutational MANOVA (through the adonis2() function) is a technique analog to ANOVA but applicable to multidimensional data, that is, when our response variable is not a single variable but an array of many variables (as is our case in which every sample is characterized by the abundance of several functions). We can test, thus, if the difference in the functional composition of the metatranscriptomes is explained by the grouping of these samples in different categories.

For example, you can test if the addition of an inoculum explains a significant portion of between-samples dissimilarity:

.. code-block:: r

    adonis2(kotab_bc~inoculant,data = metadata,permutations = 10000)

-----------------------------------------
Differential analysis
-----------------------------------------
So far we have analyzed in different ways how similar or different each pair of samples look like based on their transcript abundances for a set of metabolic functions. And we have done this by computing an index of similarity (the Bray-Curtis index) which provides a value between 0 and 1 informing on the similarity based on all functions.

Another big family of possible analyses is Differential abundance, that is, the detection of which functions are differentially abundant in a set of samples compared to another set of samples. The statistics behind these methods is out of the scope of this guide. It is sufficient at this stage to know that you can test if a specific function is differentially abundant in the samples inoculated with a bacterium in comparison to the samples without inoculum. And this can be done sistematically for all functions and detect those with significant differences.

This is precisely what you do with the following code:

.. code-block:: r

    dds<-DESeqDataSetFromMatrix(countData = t(round(kotab_matrix)),colData = metadata,design = ~ inoculant+substrate+inoculant:substrate)
    dds <- DESeq(dds)
    resultsNames(dds)
    resLFC <- lfcShrink(dds, coef="inoculant_vs_none", type="apeglm")
    res<-resLFC %>%
      as.data.frame() %>%
      arrange(desc(abs(log2FoldChange))) %>%
      rownames_to_column("KO") %>%
      left_join(distinct(select(kegg_annot,KO,Description)),by="KO") %>%
      filter(abs(log2FoldChange)>=2 & padj<=0.05)

An object called res will store the results of the statistical tests which are significant. Specifically we have tested for significant differences betweem the samples with and without inoculum, as this is the main effect detected above through PCoA and Permutational MANOVA. The main variables to consider here are:

#.    KO the KO being tested
#.    log2FoldChange the log2-transformed fold change: that is the log2-transformed ratio of the gene abundance in the samples with inoculum divided by the abundance in the samples without inoculum.
#.    padj the adjusted P-value

Thus the KOs with a log2FoldChange larger than 0 correspond to KOs overexpressed by the addition of the inoculum while the ones smaller than 0 correspond to KOs underexpressed by the addition of the inoculum.

-----------------------------------------
Visualizing the abundance of a specific KO
-----------------------------------------
With the following code you can visualize the abundance of a specific KO along the dataset. You just need to change the KO identifier in the line starting as ggplot(...). Try out one or several of the KOs annotated as ribosomal proteins:

.. code-block:: r

    toplot<-kotab_matrix %>%
      as.data.frame() %>%
      rownames_to_column("replicate_name") %>%
      left_join(metadata,by="replicate_name")

    ggplot(data=toplot,aes(x=substrate,y=K07040,color=inoculant,fill=inoculant)) +
      geom_violin(alpha=0.5,draw_quantiles = 0.5,scale = "width") +
      geom_jitter(position=position_jitterdodge()) +
      #scale_y_log10() +
      theme_minimal()

In fact the percentage of transcripts dedicated by a cell to ribosomal proteins is a good proxy of its growth state. We can compute the percentage of transcripts dedicated by the entire native community to ribosomal proteins as a proxy for its growth state:

.. code-block:: r

    rp_df<-kegg_annot %>%
      distinct(KO,Description) %>%
      filter(grepl("subunit ribosomal protein",Description) & KO %in% colnames(kotab_matrix))

    activity_df<-data.frame(total_transcripts=rowSums(kotab_matrix),
                            rp_transcripts=rowSums(kotab_matrix[,which(colnames(kotab_matrix) %in% rp_df$KO)])) %>%
      mutate(rp_perc_transcripts=100*rp_transcripts/total_transcripts) %>%
      rownames_to_column("replicate_name") %>%
      left_join(metadata,by="replicate_name")

The data produced above is stored in a data frame called activity_df. You will find a variable named rp_perc_transcripts which is the percentage of transcripts dedicated to ribosomal proteins, our proxy for growth.
