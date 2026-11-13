::: {.callout-note .callout-glance appearance="simple" collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**[The data](#the-data)**

```r
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
suppressPackageStartupMessages(library(DESeq2))
theme_set(theme_bw())

samples <- read_csv("data/samples.csv", show_col_types = FALSE) %>%
  mutate(Treatment = factor(Treatment, levels = c("PBS", "LPS")))
genes <- read_csv("data/genes.csv.gz", show_col_types = FALSE)
counts <- read_csv("data/counts.csv.gz", show_col_types = FALSE)

# the three strains with enough reads (Part 2), as a gene x sample matrix
strains <- c("I48", "YL32", "YL58")
count_matrix <- counts %>%
  filter(locus_tag %in% genes$locus_tag[genes$genome %in% strains]) %>%
  tibble::column_to_rownames("locus_tag") %>%
  as.matrix()
storage.mode(count_matrix) <- "integer"   # DESeq2 expects whole numbers
gene_strain <- genes$genome[match(rownames(count_matrix), genes$locus_tag)]
cat(nrow(count_matrix), "genes from", length(strains), "strains in", ncol(count_matrix), "samples\n")
```

**[Step 1: size factors](#step-1-size-factors)**

```r
sf <- data.frame(sample_id = colnames(count_matrix),
                 global = estimateSizeFactorsForMatrix(count_matrix))
for (s in strains) {
  sf[[s]] <- estimateSizeFactorsForMatrix(count_matrix[gene_strain == s, ])
}

sf %>%
  pivot_longer(all_of(strains), names_to = "strain", values_to = "per_strain") %>%
  left_join(samples, by = "sample_id") %>%
  ggplot(aes(global, per_strain, colour = Treatment)) +
  geom_abline(linetype = 2) +
  geom_point(size = 2) +
  facet_wrap(~ strain) +
  labs(x = "Global size factor", y = "Per-strain size factor")
```

**[Step 2: DESeq2 with per-strain size factors](#step-2-deseq2-with-per-strain-size-factors)**

```r
norm_factors <- matrix(1, nrow(count_matrix), ncol(count_matrix),
                       dimnames = dimnames(count_matrix))
for (s in strains) {
  norm_factors[gene_strain == s, ] <- matrix(sf[[s]], sum(gene_strain == s),
                                             ncol(count_matrix), byrow = TRUE)
}

coldata <- samples %>% tibble::column_to_rownames("sample_id")
dds <- DESeqDataSetFromMatrix(count_matrix, coldata[colnames(count_matrix), , drop = FALSE],
                              design = ~ Treatment)
normalizationFactors(dds) <- norm_factors
dds <- DESeq(dds, quiet = TRUE)

res <- results(dds, contrast = c("Treatment", "LPS", "PBS"), alpha = 0.01) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("locus_tag") %>%
  mutate(strain = gene_strain[match(locus_tag, rownames(count_matrix))])
res %>% arrange(padj) %>% head()
```

**[Step 3: look at the results](#step-3-look-at-the-results)**

```r
res %>%
  filter(!is.na(padj)) %>%
  mutate(call = case_when(
    padj < 0.01 & log2FoldChange > 1 ~ "up in LPS",
    padj < 0.01 & log2FoldChange < -1 ~ "down in LPS",
    TRUE ~ "not significant"
  )) %>%
  ggplot(aes(log2FoldChange, -log10(padj), colour = call)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  facet_wrap(~ strain) +
  scale_colour_manual(values = c("up in LPS" = "firebrick", "down in LPS" = "steelblue",
                                 "not significant" = "grey70")) +
  labs(x = "log2 fold change, LPS vs PBS", y = "-log10 adjusted p-value", colour = NULL)
```

:::
