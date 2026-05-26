-------------------------------------------------------------
Metatranscriptomics without Metagenomics (defined Community)
-------------------------------------------------------------
Protocol provided by Anna Sintsova.

Metatranscriptomic data arising from defined communities (i.e. community, whose composition is known) can be analysed in a way that's similar to traditional RNASeq with a few key differences. This instruction covers exploratory data analysis, differential expression, and functional enrichment.

Setting up the R Environment and downloading the Data
-----------------------------------------------------

In the first step we install the required R packages:

.. code-block:: r

    # Install all R packages required for the metatranscriptomics session.
    # Source this script once before running the course notebooks:
    #   source("install_packages.R")

    # --- CRAN packages ---
    cran_packages <- c(
      "tidyverse",
      "RColorBrewer",
      "scales",
      "ggplot2",
      "plotly",
      "pheatmap",
      "compositions",
      "data.table",
      "httr",
      "htmltools"
    )

    install.packages(cran_packages)

    # --- Bioconductor packages ---
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    bioc_packages <- c(
      "DESeq2",
      "GO.db",
      "topGO",
      "enrichplot",
      "clusterProfiler",
      "KEGGREST"
    )

    BiocManager::install(bioc_packages)

For the analysis, we will need :download:`these three datasets <../downloads/metat_defined_community_data.tar.gz>`:

=================================   ===================================================================
``METAT_raw_data.RData``                Raw counts, sample metadata, and functional annotations
``METAT_DE_results_I48GO.RData``        Differential expression results with GO annotations
``METAT_DE_results.RData``              Differential expression results (alternative, currently unused)
=================================   ===================================================================

You can download the lecture slides :download:`here <../downloads/METAT_lecture_slides.pdf>`

Part 1: Exploratory Data Analysis
---------------------------------

Goals
^^^^^
1.  Understand which strains dominate the transcriptome

2.  Check if we have adequate sequencing depth per strain

3.  Learn about different normalization strategies

4.  Use PCA and correlation heatmaps for quality control

1.1 Library composition overview
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In metatranscriptomics, understanding the compositional structure of your data is critical. Unlike single-organism RNAseq, changes in microbial community composition can confound gene expression patterns. Before we look at differential expression, we want to visualize:

1.  **Which strains contribute most to the transcriptome?**
2.  **Does treatment (LPS vs PBS) affect strain abundances?**

.. code-block:: r

    rm(list=ls())

    library(tidyverse)      # For data manipulation and ggplot2
    library(RColorBrewer)   # For color palettes
    library(scales)         # For percentage formatting
    library(ggplot2)
    theme_set(theme_bw(base_size = 12))

Load data
~~~~~~~~~

Loading raw feature counts, functional annotations, and sample metadata

.. code-block:: r

    #raw_counts <- read_csv("../data/METAT_raw_counts.csv.gz")
    #sample_metadata <- read_csv("../data/METAT_sample_metadata.csv")
    #functional_annotation <- read_table("../data/METAT_functional_annotation.tsv.gz")

    # Load raw_counts, sample_data, and annotations objects

    load("../data/METAT_raw_data.RData")

    # Preview the data structure
    cat("Count matrix dimensions:", dim(raw_counts), "\n")
    head(raw_counts[, 1:5])

Examine sample metadata
~~~~~~~~~~~~~~~~~~~~~~~
How many samples? How many conditions?

.. code-block:: r

    sample_data

Examine functional annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-   How many genes do we have per genome?

-   How many genes have a KO?

.. code-block:: r

    # Clean KO column: We convert '-', 'nan', and empty strings to a proper R 'NA' value.
    functional_annotation <- annotations %>%
      mutate(KO = case_when(
        KO == "-" ~ NA_character_,
        KO == "nan" ~ NA_character_,
        KO == "" ~ NA_character_,
        TRUE ~ KO
      ))

    functional_annotation %>% head()

    # Groupby genome, and calculate stats

    annotations_summary <- functional_annotation %>% group_by(genome) %>%
      summarise(
        # Count unique locus_tags per genome
        count_genes = n_distinct(locus_tag),

        # Count unique KOs (na.rm=TRUE automatically ignores the NAs we created)
        #count_unique_KOs =

        # Calculate percentage:
        # logic: (Count of unique genes that represent a KO / Total unique genes) * 100
        # We filter locus_tags where KO is NOT NA to get the numerator
        perc_genes_with_KO = round((n_distinct(locus_tag[!is.na(KO)]) / n_distinct(locus_tag)) * 100, 1)
      )

Examine count data
~~~~~~~~~~~~~~~~~~
What do count distributions look like?

.. code-block:: r

    raw_counts %>% head()

    # Create an interactive histogram

    raw_counts$AU647 %>% hist(breaks=100)

Which strains are actually active in our samples?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We need to aggregate gene-level counts to **strain-level counts**. This requires several steps:

-   Merge with functional_annotations to get gene-to-strain mapping

-   The count matrix is in "wide" format (genes as rows, samples as columns). We need to:

    1.  Convert to "long" format for easier aggregation: pivot the count data long
    2.  Calculate total counts per sample (library size)
    3.  Calculate the total counts per strain per sample.
    4.  Calculate each strain's percentage of the transcriptome

-   Aggregation: we use `group_by` combined with `mutate`. `mutate` adds a new column to the existing data frame while respecting the grouping. This calculates the totals while keeping the structure intact

.. code-block:: r
    # Join counts with genome IDs
    raw_counts_with_strain <- raw_counts %>%
      left_join(functional_annotation %>% select(locus_tag, genome), by = "locus_tag")

    # Transform from Wide (samples as columns) to Long (one row per observation)
    raw_counts_long <- raw_counts_with_strain %>%
      pivot_longer(
        cols = -c(locus_tag, genome), # Select all columns EXCEPT ID and genome to pivot
        names_to = "sample_id",
        values_to = "fcnts"
      )

    raw_counts_long %>% head()

.. code-block:: r
    raw_counts_processed <- raw_counts_long %>%
      # 1. Calculate Sample Total (Library Size)
      group_by(sample_id) %>%
      mutate(sample_total = sum(fcnts)) %>%

      # 2. Calculate Genome Total (Sum of counts for specific genome in specific sample)
      group_by(sample_id, genome) %>%
      mutate(genome_total = sum(fcnts))  %>% ungroup() %>%


      # 3. Calculate Raw Percentage
      mutate(genome_perc = round(genome_total / sample_total, 4))

    raw_counts_processed %>% head()

-   Visualising the composition: we extract the genome-level statistics and join with the sample metadata (`sample_data`).

.. code-block:: r

    # Create a summarized table of just genome stats (removing individual gene IDs)
    genome_counts <- raw_counts_processed %>%
      select(sample_id, genome, sample_total, genome_total, genome_perc) %>%
      distinct()

    # Join metadata and filter
    oligo_metat_composition<- genome_counts %>%
      left_join(sample_data, by = "sample_id") %>% arrange(Treatment)
    oligo_metat_composition

-   Visualize transcriptome composition. **Important**: the transcriptome composition does not directly represent taxonomic composition. Why?

.. code-block:: r
    new_order <- sample_data %>% arrange(Treatment)%>% pull(sample_id)
    oligo_metat_composition$sample_id <- factor(
      oligo_metat_composition$sample_id,
      levels = new_order)


    # Define custom colors
    syncom_colors <- c(
      "YL32" = "#149AB3",
      "KB18" = "#616161",
      "I48"  = "#C26215",
      "YL27" = "#F59C46",
      "YL45" = "#E30C4B",
      "I46"  = "#2C5D52",
      "YL58" = "#163A1A",
      "YL31" = "#099334",
      "YL2"  = "#282E68",
      "I49"  = "#3DB077",
      "YL44" = "#AC7FB6",
      "KB1"  = "#42AB34"
    )

    p <- ggplot(oligo_metat_composition, aes(x = sample_id, y = genome_perc, fill = genome)) +
      geom_col(position = "stack") +

      # Apply custom colors
      scale_fill_manual(values = syncom_colors) +

      # Add the vertical dashed line (x = 5.5)
      geom_vline(xintercept = 6.5, linetype = "dashed", color = "black", linewidth = 0.8) +

      # Add Annotations (LPS and PBS)
      # Note: Check x-coordinates if sample_id is categorical.
      # If categorical, 3 and 8 correspond to the 3rd and 8th items on the axis.
      annotate("text", x = 3, y = 1.03, label = "LPS", fontface = "bold") +
      annotate("text", x = 8, y = 1.03, label = "PBS", fontface = "bold") +

      # Formatting
      theme_minimal() +
      labs(y = "Transcriptome abundance", x = NULL) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank()
      )

    # Display static plot
    p

    # Convert to interactive Plotly (optional)
    ggplotly(p)

To estimate taxonomic composition, we could normalize by genome length, or (better) use other tools /sequencing methods (i.e. 16S, mOTUs on transcriptomic data).

1.2 Sequencing Depth and Sample Quality Assessment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In metatranscriptomics, total library size isn't enough - we need adequate coverage **per strain** to detect differential expression. A sample with 10M total reads might have:

-   8M reads from one abundant strain (excellent coverage)

-   100K reads from a rare strain (poor coverage)

For today, we want **≥500K reads per strain** for reliable DE analysis

.. code-block:: r

    # Define threshold for adequate coverage
    coverage_threshold <- 5e5  # 500,000 reads

    strain_depth <- oligo_metat_composition %>%
      mutate(
        adequate_coverage = genome_total >= coverage_threshold,
        counts_millions = genome_total / 1e6
      )

    strain_depth

Samples with adequate coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-   This table shows how many replicates per treatment have sufficient coverage for each strain.

.. code-block:: r

    coverage_summary <- strain_depth %>%
      group_by(genome, Treatment) %>%
      summarize(
        n_samples = n(),
        n_adequate = sum(adequate_coverage),
        mean_counts = round(mean(genome_total),2),
        min_counts = min(genome_total),
        max_counts = max(genome_total),
        .groups = "drop"
      ) %>%
      mutate(
        mean_millions = mean_counts / 1e6,
        min_millions = min_counts / 1e6,
        max_millions = max_counts / 1e6
      ) %>%


      select(genome, Treatment, n_adequate, mean_millions) %>%
      pivot_wider(
        names_from = Treatment,
        values_from = c(n_adequate, mean_millions),
        names_sep = "_"
      ) %>%
      arrange(desc(mean_millions_LPS))

    coverage_summary

    coverage_bubble_data <- strain_depth %>%
      group_by(genome, Treatment) %>%
      summarize(
        n_adequate = sum(adequate_coverage),
        n_total = n(),
        .groups = "drop"
      )

    ggplot(coverage_bubble_data,
           aes(x = Treatment, y = genome, size = n_adequate, color = genome)) +
      geom_point(alpha = 0.7) +

      # Add text labels with the number
      geom_text(
        aes(label = n_adequate),
        color = "white",
        fontface = "bold",
        size = 4,
        show.legend = FALSE
      ) +

      scale_size_continuous(
        name = "n replicates\n≥500K reads",
        range = c(1, 10),
        breaks = c(0, 1, 2, 3, 4, 5, 6),
        limits = c(0, NA)
      ) +

      scale_color_manual(values = syncom_colors, guide = "none") +

      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.border = element_rect(fill = NA, color = "grey70")
      ) +

      labs(
        x = NULL,
        y = NULL,
        title = "Samples with adequate sequencing depth per strain",
        subtitle = "Number shown = replicates ≥500K reads"
      )

Data Normalization for Exploratory Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Why normalize?
~~~~~~~~~~~~~~

Before we perform PCA and other exploratory analyses, we need to normalize our count data to account for:

1.  **Gene length bias:** Longer genes generate more reads even at equal expression
2.  **Library size differences:** Samples with more total reads have higher counts
3.  **Compositionality:** In (meta)transcriptomics, counts are relative (can't measure "absolute" expression)

**Important distinction:** - **For visualization/PCA:** We'll use TPM → log2 or CLR transformation - **For differential expression:** DESeq2 does its own normalization (median-of-ratios with taxon-specific scaling)

Method 1: TPM (Transcripts Per Million)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**TPM normalization:** 1. Divide counts by gene length (in kb) → reads per kilobase (RPK) 2. Sum all RPK values in the sample → "per million" scaling factor 3. Divide each RPK by this scaling factor and multiply by 10\^6

**Result:** TPM values are comparable: - Within samples (gene-to-gene) - Between samples (same gene across conditions) - Sum to \~1 million per sample

.. code-block:: r

    # Function to calculate TPM
    calculate_tpm <- function(counts, lengths) {
      # counts: numeric vector of read counts
      # lengths: numeric vector of gene lengths in bp

      # Step 1: Calculate RPK (reads per kilobase)
      rpk <- counts / (lengths / 1000)

      # Step 2: Calculate scaling factor (sum of all RPK / 1 million)
      scaling_factor <- sum(rpk, na.rm = TRUE) / 1e6

      # Step 3: Calculate TPM
      tpm <- rpk / scaling_factor

      return(tpm)
    }

.. code-block:: r

    # Check if all genes in count matrix have length info
    functional_annotation <- functional_annotation %>% mutate(gene_length = End - Start)
    # Add gene lengths to count data
    counts_with_length <- raw_counts_long %>%
      left_join(functional_annotation %>% select(locus_tag, gene_length),
                by = 'locus_tag')

    tpm_data <- counts_with_length %>%
      group_by(genome, sample_id) %>%
      mutate(tpm = calculate_tpm(fcnts, gene_length)) %>%
      ungroup()



    tpm_sums <- tpm_data %>%
      group_by(genome, sample_id) %>%
      summarize(total_tpm = sum(tpm, na.rm = TRUE))

    tpm_sums

Transformation Option 1: Log2 transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Log transformation is standard for:

\- Variance stabilization

\- Making data more normally distributed

**Key decision:** What pseudocount to add? (Can't take log of zero)

.. code-block:: r

    # Add pseudocount and log2 transform
    # Common choices: 0.5, 1, or 0.25
    pseudocount <- 0.5

    # Apply log2 transformation to TPM values
    tpm_data <- tpm_data %>%
      mutate(log2_tpm = log2(tpm + pseudocount))

    # Check distribution across samples
    ggplot(tpm_data, aes(x = log2_tpm, color = sample_id)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(
        title = "Distribution of log2(TPM + 0.5) values",
        x = "log2(TPM + 0.5)",
        y = "Density"
      )

Transformation Option 2: CLR (Centered Log-Ratio)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CLR transformation is **compositionally-aware** - important for (meta)transcriptomics!

**Why CLR?**

\- Explicitly accounts for compositional nature of sequencing data

\- Preserves relative relationships

\- Better for PCA of compositional data (?)

**How it works:**

1\. Replace zeros with a small value (geometric mean cannot handle zeros)

2\. Calculate geometric mean of all genes in a sample

3\. Take log-ratio of each gene to the geometric mean

.. code-block:: r

    # Function to calculate CLR transformation
    calculate_clr <- function(x, pseudocount = 0.5) {
      # x: vector of TPM values for one sample

      # Replace zeros with small pseudocount
      x_no_zero <- x + pseudocount

      # Calculate geometric mean
      gm <- exp(mean(log(x_no_zero)))

      # CLR = log(value / geometric_mean)
      clr <- log(x_no_zero / gm)

      return(clr)
    }

    # Apply CLR transformation per sample
    tpm_data <- tpm_data %>%
      group_by(genome, sample_id) %>%
      mutate(clr_tpm = as.numeric(compositions::clr(tpm + 0.5)),
             clr_fcnt = calculate_clr(fcnts, pseudocount=0.5),
             clr_custom = calculate_clr(tpm, pseudocount = 0.5)) %>%  # pseudocount for zeros
      ungroup()

    tpm_data %>% head()

    # Filter to one genome and one sample, then reshape for comparison
    comparison_data <- tpm_data %>%
      filter(genome == "YL58", sample_id == "AU658") %>%
      select(locus_tag, tpm, log2_tpm,  clr_custom, clr_fcnt) %>%
      pivot_longer(cols = c(log2_tpm, tpm, clr_fcnt),
                   names_to = "normalization",
                   values_to = "value")

    # Create comparison plot
    ggplot(comparison_data, aes(x = value)) +
      geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
      facet_wrap(~normalization, scales = "free") +
      theme_bw() +
      labs(title = "Comparison of Normalization Methods",
           subtitle = "Genome_A, Sample_1",
           x = "Normalized Value",
           y = "Count")

    genome = 'YL58'
    # Scatter plot comparing two normalizations
    tpm_data %>%
      filter(genome == 'I48', sample_id == 'AU649') %>%
      ggplot(aes(x = log2_tpm, y = clr_fcnt
                 )) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      theme_bw() +
      labs(title = "Comparison: log2(TPM) vs CLR(TPM)",
           subtitle = "Genome_A, Sample_1",
           x = "log2(TPM + 1)",
           y = "CLR(TPM)")

Principal Component Analysis (PCA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PCA is a dimensionality reduction technique that helps visualize patterns in high-dimensional gene expression data. It transforms the data into principal components (PCs) that capture the most variation, allowing us to see how samples cluster based on their overall expression profiles.

We'll perform PCA on CLR-normalized data for a single genome and color samples by a metadata variable to identify biological patterns.

The percentage values on each axis indicate how much of the total variance is explained by that component.

.. code-block:: r

    # Prepare data: pivot to wide format (samples as rows, genes as columns)
    pca_matrix <- tpm_data %>%
      filter(genome == 'I48') %>%
      select(sample_id, locus_tag, clr_tpm) %>%
      pivot_wider(names_from = locus_tag, values_from = clr_tpm) %>%
      column_to_rownames("sample_id")

    # Run PCA
    pca_result <- prcomp(pca_matrix, center = TRUE, scale. = FALSE)

    # Extract PC scores and merge with metadata
    pca_scores <- as.data.frame(pca_result$x) %>%
      rownames_to_column("sample_id") %>%
      left_join(sample_data, by = "sample_id")

    # Plot
    ggplot(pca_scores, aes(x = PC1, y = PC2, color = Treatment)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(title = "PCA of Gene Expression (CLR normalized)",
           subtitle = 'I48',
           x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
           y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))

    # Prepare data: pivot to wide format (samples as rows, genes as columns)
    pca_matrix <- tpm_data %>%
      filter(genome == 'I48') %>%
      select(sample_id, locus_tag, clr_fcnt) %>%
      pivot_wider(names_from = locus_tag, values_from = clr_fcnt) %>%
      column_to_rownames("sample_id")

    # Run PCA
    pca_result <- prcomp(pca_matrix, center = TRUE, scale. = FALSE)

    # Extract PC scores and merge with metadata
    pca_scores <- as.data.frame(pca_result$x) %>%
      rownames_to_column("sample_id") %>%
      left_join(sample_data, by = "sample_id")

    # Plot
    ggplot(pca_scores, aes(x = PC1, y = PC2, color = Treatment)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(title = "PCA of Gene Expression (CLR normalized)",
           subtitle = 'I48',
           x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
           y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))

    run_pca <- function(data, sample_data, genome_name, norm_column, color_by, top_n_genes = NULL) {

      # Prepare data: filter and pivot to wide format
      data_filtered <- data[data$genome == genome_name, ]
      data_selected <- data_filtered[, c("sample_id", "locus_tag", norm_column)]
      pca_wide <- pivot_wider(data_selected, names_from = locus_tag, values_from = all_of(norm_column))
      print(dim(pca_wide))
      pca_matrix <- as.data.frame(pca_wide)
      rownames(pca_matrix) <- pca_matrix$sample_id
      pca_matrix$sample_id <- NULL

      # Filter for most variable genes if specified
      if (!is.null(top_n_genes)) {
        gene_variance <- apply(pca_matrix, 2, var)
        top_genes <- names(sort(gene_variance, decreasing = TRUE)[1:top_n_genes])
        pca_matrix <- pca_matrix[, top_genes]
      }

      # Run PCA
      pca_result <- prcomp(pca_matrix)

      # Extract PC scores and merge with metadata
      pca_scores <- as.data.frame(pca_result$x)
      pca_scores$sample_id <- rownames(pca_scores)
      pca_scores <- merge(pca_scores, sample_data, by = "sample_id")

      # Plot
      ggplot(pca_scores, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
        geom_point(size = 3) +
        theme_bw() +
        labs(title = paste0("PCA of Gene Expression (", norm_column, ")"),
             subtitle = paste0(genome_name,
                              ifelse(!is.null(top_n_genes),
                                    paste0(" - Top ", top_n_genes, " variable genes"), "")),
             x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
             y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))
    }

**Usage examples:**

.. code-block:: r

    # PCA with all genes
    run_pca(tpm_data, sample_data, "I48", "log2_tpm", "Treatment", top_n_genes=500)

    # PCA with top 500 most variable genes
    run_pca(tpm_data, sample_data, "I48", "clr_tpm", "Treatment", top_n_genes = 500)

    # Different genome and normalization
    run_pca(tpm_data, sample_data, "I48", "log2_tpm", "Treatment", top_n_genes = 3000)

Filtering for the most variable genes can improve visualization by reducing noise from genes with minimal expression variation.

Sample Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^

Sample-to-sample correlation heatmaps help us assess the quality and consistency of our data. High correlations between biological replicates indicate good technical reproducibility, while outlier samples with low correlations may indicate problems during library preparation or sequencing. This visualization is useful for quality control before downstream analyses.

.. code-block:: r

    library(pheatmap)


    # Prepare correlation matrix for a single genome
    correlation_matrix <- tpm_data %>%
      filter(genome == "I48") %>%
      select(sample_id, locus_tag, log2_tpm) %>%
      pivot_wider(names_from = sample_id, values_from = log2_tpm) %>%
      column_to_rownames("locus_tag") %>%
      cor(method = "spearman")

    # Optional: prepare annotation for metadata
    annotation_col <- sample_data %>%
      select(sample_id, Treatment) %>%  # adjust columns as needed
      column_to_rownames("sample_id")

    # Create heatmap
    pheatmap(correlation_matrix,
             annotation_col = annotation_col,
             annotation_row = annotation_col,
             color = colorRampPalette(c("white", "darkblue"))(50),,
             breaks = seq(0.7, 1, length.out = 51),  # adjust range as needed
             main = "Sample-to-Sample Correlation",
             display_numbers = TRUE,
             clustering_method = 'ward.D2',
             number_format = "%.2f",
             fontsize_number = 8)

- Investigate different clustering methods - how do they affect the grouping?
- Investigate whether the patterns are similar for different strains in the community?


Differential Expression
-----------------------

Differential expression (DE) analysis in metatranscriptomics identifies
genes that are significantly differentially expressed between
conditions. Unlike single-organism transcriptomics, metatranscriptomic
data requires special consideration:

**Key challenge:** Different organisms have different abundances, which
can confound normalization. If one organism blooms in a condition, its
genes will appear highly expressed simply due to increased biomass, not
true differential regulation.

**Solution:** We normalize **within each genome** separately, then
perform DE analysis on the combined dataset. This approach:

1.  Calculates size factors for each genome independently
2.  Normalizes counts within each genome
3.  Runs differential expression on the combined normalized data

This workflow is adapted from *Hierarchical normalization for
differential expression analysis in metatranscriptomics* (Jorth et al.).

Setup:

.. code-block:: r

    rm(list=ls())
    library(RColorBrewer)   # For color palettes
    library(scales)         # For percentage formatting
    library(plotly)         # For interactive plots (optional)

    library(DESeq2)
    library(tidyverse)
    library(data.table)

    theme_set(theme_bw(base_size = 12))


**Expected data structure:** - `count_data`: columns for locus_tag/ID,
genome, and one column per sample with raw counts - `sample_data`:
sample_id, condition/treatment variables (e.g., Treatment) -
`annotations`: functional annotations for genes (optional, for
interpretation)

Helper Functions
^^^^^^^^^^^^^^^^

These functions implement the within-genome normalization approach:

.. code-block:: r

    # Align sample data and count data
    align_samples_and_counts <- function(count_data, sample_data, conditions,
                                         gene_col, sample_col) {
      count_data <- count_data %>% as.data.frame()
      rownames(count_data) <- count_data[[gene_col]]

      # Filter samples present in both datasets
      sample_data <- sample_data %>% filter(sample_id %in% colnames(count_data))

      # Create combined group variable from conditions
      conditions <- unlist(strsplit(conditions, ","))
      sample_data <- sample_data %>% unite("group", all_of(conditions), remove = FALSE)

      # Reorder count data to match sample data
      count_data <- count_data[, sample_data[[sample_col]]]

      return(list("count_data" = count_data, "sample_data" = sample_data))
    }

    # Normalize counts for a single genome using DESeq2
    normalize_single_genome <- function(raw_data) {
      count_data <- raw_data$count_data
      col_data <- raw_data$sample_data
      count_data_row <- nrow(count_data)
      count_data <- round(count_data)  # Ensure integer counts

      # Check if there are non-zero rows
      if (sum(rowSums(count_data == 0) == 0) != 0) {

        # Create DESeq2 object
        dds <- DESeqDataSetFromMatrix(countData = count_data,
                                       colData = col_data,
                                       design = ~group)
        colData(dds)$group <- factor(colData(dds)$group,
                                      levels = unique(colData(dds)$group))

        # Filter low-count genes
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep, ]

        # Calculate size factors
        dds <- DESeq(dds, quiet = TRUE)

        # Normalize counts using size factors
        norm_count <- count_data / rep(dds@colData@listData$sizeFactor,
                                       each = count_data_row)
        return(norm_count)
      } else {
        return(data.frame())
      }
    }

    # Normalize within each genome
    normalize_by_genome <- function(count_data, sample_data, conditions,
                                   gene_col, sample_id_col, genome_col) {
      genomes <- count_data[[genome_col]] %>% unique()
      cat("Normalizing within each genome. Number of genomes:", length(genomes), "\n")

      norm_list <- list()
      for (i in seq_along(genomes)) {
        cat("Processing:", genomes[[i]], "\n")

        # Extract data for this genome
        genome_data <- count_data %>%
          filter(!!sym(genome_col) == genomes[[i]]) %>%
          select(-!!sym(genome_col))

        # Align and normalize
        raw_data <- align_samples_and_counts(genome_data, sample_data,
                                             conditions, gene_col, sample_id_col)
        norm_list[[i]] <- normalize_single_genome(raw_data) %>%
          rownames_to_column(gene_col)
      }

      cat("Normalization complete!\n")
      return(rbindlist(norm_list, fill = TRUE))
    }

    # Run DE analysis on normalized data
    run_de_analysis <- function(norm_counts, sample_data, conditions,
                               gene_col, sample_id_col) {

      # Prepare normalized counts
      rownames(norm_counts) <- norm_counts[[gene_col]]
      norm_counts <- norm_counts %>%
        select(-!!sym(gene_col)) %>%
        round() %>%
        rownames_to_column(gene_col)

      # Align samples
      raw_data <- align_samples_and_counts(norm_counts, sample_data,
                                           conditions, gene_col, sample_id_col)

      # Create DESeq2 object with pre-normalized data
      dds <- DESeqDataSetFromMatrix(countData = raw_data$count_data,
                                     colData = raw_data$sample_data,
                                     design = ~group)
      colData(dds)$group <- factor(colData(dds)$group,
                                   levels = unique(colData(dds)$group))

      # Set normalization factors to 1 (data already normalized)
      normFactors <- matrix(1, ncol = ncol(raw_data$count_data),
                            nrow = nrow(raw_data$count_data))
      normalizationFactors(dds) <- normFactors

      cat("Running differential expression analysis...\n")
      dds <- DESeq(dds)
      return(dds)
    }

    # Main wrapper function
    run_metatranscriptomics_de <- function(count_data, sample_data, conditions,
                                           gene_col, sample_id_col, genome_col) {
      # Step 1: Normalize by genome
      norm_counts <- normalize_by_genome(count_data, sample_data, conditions,
                                         gene_col, sample_id_col, genome_col)

      # Step 2: Run DESeq2
      dds <- run_de_analysis(norm_counts, sample_data, conditions,
                             gene_col, sample_id_col)

      return(list("norm_counts" = norm_counts, "dds" = dds))
    }

Define Analysis Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: r

    # Define conditions to compare
    conditions <- "Treatment"  # Can be multiple: "Treatment,Timepoint"

    # Statistical thresholds
    lfc_threshold <- 0.0    # Log2 fold-change threshold
    alpha <- 0.01           # Adjusted p-value threshold

    # Column names in your data
    gene_col <- "locus_tag"           # Gene identifier column
    sample_col <- "sample_id"  # Sample identifier column
    genome_col <- "genome"     # Genome/taxon identifier column

Filter Data (Optional)
^^^^^^^^^^^^^^^^^^^^^^

If you want to analyze specific genomes or samples:

.. code-block:: r

    # Example: Select specific genomes
    genomes_to_analyze <- c("YL32", "I48", "YL58")

    count_data_filtered <- count_data %>%
      filter(genome %in% genomes_to_analyze)

    # Use all samples or filter by condition
    sample_data_filtered <- sample_data  # Use all samples

Run Analysis
^^^^^^^^^^^^

.. code-block:: r

    # Run the complete workflow
    de_results <- run_metatranscriptomics_de(
      count_data = count_data_filtered,
      sample_data = sample_data_filtered,
      conditions = conditions,
      gene_col = gene_col,
      sample_id_col = sample_col,
      genome_col = genome_col
    )

    # Extract results
    dds <- de_results$dds
    norm_counts <- de_results$norm_counts

Extract Results
^^^^^^^^^^^^^^^

Get Pairwise Comparisons
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: r

    # Define contrasts: compare groups containing these patterns
    pattern1 <- "LPS"   # Treatment group
    pattern2 <- "PBS"   # Control group

    # Identify all groups matching patterns
    sample_groups <- colData(dds)
    first_contrasts <- unique(sample_groups$group[grepl(pattern1, sample_groups$group)])
    second_contrasts <- unique(sample_groups$group[grepl(pattern2, sample_groups$group)])
    first_contrasts <- setdiff(first_contrasts, second_contrasts)

    # Create all pairwise comparisons
    contrasts_to_test <- expand.grid(
      treatment = first_contrasts,
      control = second_contrasts
    )

    print(contrasts_to_test)

Extract DE Results
~~~~~~~~~~~~~~~~~~

.. code-block:: r
    # Function to extract results for one comparison
    get_de_results <- function(dds, condition1, condition2, lfc_thresh, alpha_thresh) {
      res <- results(dds,
                     contrast = c("group", condition1, condition2),
                     alpha = alpha_thresh,
                     lfcThreshold = lfc_thresh)

      res_df <- as.data.frame(res) %>%
        rownames_to_column("ID") %>%
        mutate(contrast = paste0(condition1, "_vs_", condition2)) %>%
        arrange(padj)

      return(res_df)
    }

    # Get results for all comparisons
    all_results <- list()
    for (i in 1:nrow(contrasts_to_test)) {
      c1 <- as.character(contrasts_to_test$treatment[i])
      c2 <- as.character(contrasts_to_test$control[i])

      cat("Extracting results for:", c1, "vs", c2, "\n")

      all_results[[paste0(c1, "_vs_", c2)]] <- get_de_results(
        dds, c1, c2, lfc_threshold, alpha
      )
    }

    # View first comparison
    head(all_results[[1]])

Summarize Results
~~~~~~~~~~~~~~~~~

.. code-block:: r

    # Count significant genes per comparison
    summary_df <- data.frame(
      comparison = names(all_results),
      total_genes = sapply(all_results, nrow),
      significant = sapply(all_results, function(x) sum(x$padj < alpha, na.rm = TRUE)),
      upregulated = sapply(all_results, function(x) sum(x$padj < alpha & x$log2FoldChange > 0, na.rm = TRUE)),
      downregulated = sapply(all_results, function(x) sum(x$padj < alpha & x$log2FoldChange < 0, na.rm = TRUE))
    )

    print(summary_df)

Visualization
~~~~~~~~~~~~~

.. code-block:: r

    # Function to create volcano plot
    plot_volcano <- function(results_df, comparison_name, alpha_thresh = 0.01, lfc_thresh = 1) {

      results_df <- results_df %>%
        mutate(
          significance = case_when(
            padj < alpha_thresh & log2FoldChange > lfc_thresh ~ "Upregulated",
            padj < alpha_thresh & log2FoldChange < -lfc_thresh ~ "Downregulated",
            TRUE ~ "Not significant"
          )
        )

      ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = c("Upregulated" = "red",
                                       "Downregulated" = "blue",
                                       "Not significant" = "gray")) +
        geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(alpha_thresh), linetype = "dashed", color = "black") +
        theme_bw() +
        labs(title = paste("Volcano Plot:", comparison_name),
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value") +
        theme(legend.position = "bottom")
    }

    # Plot first comparison
    plot_volcano(all_results[[1]], names(all_results)[1], alpha_thresh = alpha, lfc_thresh = 1)

**MA Plot**

.. code-block:: r

    # MA plot for first comparison
    plotMA(results(dds,
                   contrast = c("group",
                              as.character(contrasts_to_test$treatment[1]),
                              as.character(contrasts_to_test$control[1])),
                   alpha = alpha),
           main = paste("MA Plot:", names(all_results)[1]),
           ylim = c(-5, 5))

**Top DE Genes Heatmap**

.. code-block:: r

    library(pheatmap)

    # Get top 50 DE genes from first comparison
    top_genes <- all_results[[1]] %>%
      filter(padj < alpha) %>%
      arrange(padj) %>%
      head(50) %>%
      pull(ID)

    # Extract normalized counts for these genes
    if (length(top_genes) > 0) {
      heatmap_data <- norm_counts %>%
        filter(locus_tag %in% top_genes) %>%
        column_to_rownames("locus_tag")

      # Create annotation
      annotation_col <- sample_data_filtered %>%
        select(sample_id, Treatment) %>%
        column_to_rownames("sample_id")

      # Plot heatmap
      pheatmap(heatmap_data,
               annotation_col = annotation_col,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               main = "Top 50 Differentially Expressed Genes",
               fontsize_row = 6,
               show_rownames = TRUE)
    }

Export Results
~~~~~~~~~~~~~~

.. code-block:: r

    # Create output directory
    output_dir <- "DE_results"
    dir.create(output_dir, showWarnings = FALSE)

    # Save normalized counts
    write.csv(norm_counts,
              file.path(output_dir, "normalized_counts.csv"),
              row.names = FALSE)

    # Save DE results for each comparison
    for (comparison in names(all_results)) {
      write.csv(all_results[[comparison]],
                file.path(output_dir, paste0(comparison, "_DESeq2_results.csv")),
                row.names = FALSE)
    }

    # Save summary
    write.csv(summary_df,
              file.path(output_dir, "DE_summary.csv"),
              row.names = FALSE)

    cat("Results saved to:", output_dir, "\n")

Key Takeaways
~~~~~~~~~~~~~

1.  **Within-genome normalization** is essential for metatranscriptomic
    data to account for differential organism abundances

2.  **DESeq2** calculates size factors within each genome, then performs
    DE analysis on the combined data

3.  **Adjusted p-values (padj)** control for multiple testing - use
    these for determining significance

4.  **Log2 fold change** indicates magnitude of differential expression
    (e.g., log2FC = 2 means 4-fold upregulation)

5.  You can visualize your results (volcano plots, MA plots, heatmaps) to
    understand global patterns