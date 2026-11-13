::: {.callout-note .callout-glance appearance="simple" collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**[The data](#the-data)**

```r
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())

samples <- read_csv("data/samples.csv", show_col_types = FALSE)
genes <- read_csv("data/genes.csv.gz", show_col_types = FALSE)
counts <- read_csv("data/counts.csv.gz", show_col_types = FALSE)
kegg_ko <- read_csv("data/kegg_ko.csv", show_col_types = FALSE)

# one row per gene per sample, with the gene's strain and length
counts_long <- counts %>%
  pivot_longer(-locus_tag, names_to = "sample_id", values_to = "count") %>%
  left_join(select(genes, locus_tag, genome, length), by = "locus_tag") %>%
  left_join(samples, by = "sample_id")

strain_colours <- c(
  YL32 = "#149AB3", KB18 = "#616161", I48 = "#C26215", YL27 = "#F59C46",
  YL45 = "#E30C4B", I46 = "#2C5D52", YL58 = "#163A1A", YL31 = "#099334",
  YL2 = "#282E68", I49 = "#3DB077", YL44 = "#AC7FB6", KB1 = "#42AB34"
)
cat("Loaded", nrow(genes), "genes from", n_distinct(genes$genome), "strains and",
    nrow(samples), "samples.\n")
cat("Ready: packages and data are loaded.\n")
```

**[Step 1: who dominates the transcriptome?](#step-1-who-dominates-the-transcriptome)**

```r
strain_inserts <- counts_long %>%
  group_by(sample_id, Treatment, genome) %>%
  summarise(inserts = sum(count), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(share = inserts / sum(inserts)) %>%
  ungroup()

ggplot(strain_inserts, aes(sample_id, share, fill = genome)) +
  geom_col() +
  facet_grid(~ Treatment, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = strain_colours) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Share of transcripts", fill = "Strain") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
```

**[Step 2: is there enough data for each strain?](#step-2-is-there-enough-data-for-each-strain)**

```r
depth <- strain_inserts %>%
  group_by(genome, Treatment) %>%
  summarise(n_adequate = sum(inserts >= 5e5), .groups = "drop")

ggplot(depth, aes(Treatment, genome, fill = n_adequate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_adequate)) +
  scale_fill_gradient(low = "grey95", high = "steelblue", limits = c(0, 6)) +
  labs(x = NULL, y = NULL, fill = "Samples with\n≥ 500k inserts")
```

**[Step 3: normalising for exploration](#step-3-normalising-for-exploration)**

```r
calculate_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  rpk / sum(rpk) * 1e6
}

expr <- counts_long %>%
  group_by(genome, sample_id) %>%
  mutate(tpm = calculate_tpm(count, length)) %>%
  ungroup()

# Each strain's TPM adds up to a million in every sample
expr %>%
  group_by(genome, sample_id) %>%
  summarise(total_tpm = sum(tpm), .groups = "drop") %>%
  head()

expr <- expr %>%
  group_by(genome, sample_id) %>%
  mutate(
    log2_tpm = log2(tpm + 0.5),
    clr = log(tpm + 0.5) - mean(log(tpm + 0.5))
  ) %>%
  ungroup()

expr %>%
  filter(genome == "I48", sample_id == "AU647") %>%
  pivot_longer(c(tpm, log2_tpm, clr), names_to = "scale") %>%
  mutate(scale = factor(scale, levels = c("tpm", "log2_tpm", "clr"))) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 60, fill = "steelblue") +
  facet_wrap(~ scale, scales = "free") +
  labs(title = "I48 in one sample: the same genes on three scales", x = NULL, y = "Genes")
```

**[Step 4: do replicates agree?](#step-4-do-replicates-agree)**

```r
strain <- "I48"

# genes with at least 10 inserts in at least 6 samples
detected <- counts_long %>%
  filter(genome == strain) %>%
  group_by(locus_tag) %>%
  summarise(n_samples = sum(count >= 10)) %>%
  filter(n_samples >= 6) %>%
  pull(locus_tag)

clr_matrix <- expr %>%
  filter(genome == strain, locus_tag %in% detected) %>%
  select(sample_id, locus_tag, clr) %>%
  pivot_wider(names_from = locus_tag, values_from = clr) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix()

# the 500 genes whose CLR values vary most between samples
top_genes <- names(sort(apply(clr_matrix, 2, var), decreasing = TRUE))[1:500]

cat(length(detected), "of", n_distinct(expr$locus_tag[expr$genome == strain]),
    strain, "genes are detected; the PCA uses the 500 most variable.\n")

pca <- prcomp(clr_matrix[, top_genes])
explained <- round(100 * summary(pca)$importance[2, 1:2])

as.data.frame(pca$x) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(samples, by = "sample_id") %>%
  ggplot(aes(PC1, PC2, colour = Treatment, label = sample_id)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(title = paste("PCA of", strain, "gene expression (CLR)"),
       x = paste0("PC1 (", explained[1], "%)"), y = paste0("PC2 (", explained[2], "%)"))
```

**[Step 5: how fast are the strains growing?](#step-5-how-fast-are-the-strains-growing)**

```r
ribosomal_kos <- kegg_ko %>%
  filter(grepl("subunit ribosomal protein", name)) %>%
  pull(KO)

# some genes have several KOs, separated by commas
genes <- genes %>%
  mutate(is_ribosomal = vapply(strsplit(coalesce(KO, ""), ","),
                               function(k) any(k %in% ribosomal_kos), logical(1)))

counts_ribo <- counts_long %>%
  left_join(select(genes, locus_tag, is_ribosomal), by = "locus_tag")

cat(length(ribosomal_kos), "ribosomal protein KOs;",
    sum(genes$is_ribosomal), "ribosomal protein genes in the community\n")

ribo <- counts_ribo %>%
  filter(genome %in% c("I48", "YL32", "YL58")) %>%
  group_by(genome, sample_id, Treatment) %>%
  summarise(ribosomal_share = sum(count[is_ribosomal]) / sum(count), .groups = "drop")

ggplot(ribo, aes(Treatment, ribosomal_share, colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  facet_wrap(~ genome, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Transcripts from ribosomal proteins", x = NULL)
```

:::
