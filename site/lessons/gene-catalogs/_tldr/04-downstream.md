::: {.callout-note .callout-glance appearance="simple" collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**[The data](#the-data)**

```r
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())

samples <- read_csv("data/tara_samples.csv", show_col_types = FALSE)
kos <- read_csv("data/tara_kos.csv", show_col_types = FALSE)
genes <- read_csv("data/tara_gene_slice.csv.gz", show_col_types = FALSE)
cat("Loaded", nrow(samples), "samples and", nrow(genes), "gene counts.\n")
```

**[Step 1: gene length](#step-1-gene-length)**

```r
ko_sums <- read_csv("data/tara_ko_sums.csv", show_col_types = FALSE)

ko_sums %>%
  filter(KO %in% c("K03040", "K03043"), grepl("_METAG$", sample)) %>%
  select(sample, KO, inserts, inserts_per_bp) %>%
  pivot_wider(names_from = KO, values_from = c(inserts, inserts_per_bp)) %>%
  summarise(
    raw_ratio = median(inserts_K03043 / inserts_K03040),
    length_normalised_ratio = median(inserts_per_bp_K03043 / inserts_per_bp_K03040)
  )
```

**[Step 2: sequencing depth and cells](#step-2-sequencing-depth-and-cells)**

```r
marker_kos <- c("K06942", "K01889", "K01887", "K01875", "K01883",
                "K01869", "K01873", "K01409", "K03106", "K03110")

mg_median <- ko_sums %>%
  filter(KO %in% marker_kos) %>%
  group_by(sample) %>%
  summarise(mg_median = median(inserts_per_bp))

mg_median %>%
  inner_join(samples, by = c("sample" = "sample_metag")) %>%
  ggplot(aes(x = sample_metag_nreads / 1e6, y = mg_median)) +
  geom_point(alpha = 0.6) +
  labs(x = "Sequencing depth (million reads)",
       y = "Median marker gene abundance\n(inserts per bp)")

per_cell <- ko_sums %>%
  inner_join(mg_median, by = "sample") %>%
  mutate(per_cell = inserts_per_bp / mg_median)

per_cell %>%
  filter(grepl("_METAG$", sample)) %>%
  group_by(KO) %>%
  summarise(median_copies_per_cell = median(per_cell)) %>%
  left_join(kos, by = "KO")
```

**[Step 4: comparing samples](#step-4-comparing-samples)**

```r
top <- read_csv("data/tara_top_genes.csv.gz", show_col_types = FALSE)

# genes per cell for each marker gene in each metagenome, as a sample x gene matrix
profile <- top %>%
  inner_join(mg_median, by = "sample") %>%
  mutate(per_cell = inserts / length / mg_median) %>%
  select(sample, gene, per_cell) %>%
  pivot_wider(names_from = gene, values_from = per_cell, values_fill = 0)

profile_matrix <- as.matrix(profile[, -1])
rownames(profile_matrix) <- profile$sample
dim(profile_matrix)

bray_curtis <- function(x) {
  n <- nrow(x)
  d <- matrix(0, n, n, dimnames = list(rownames(x), rownames(x)))
  for (i in seq_len(n)) {
    diff <- abs(sweep(x, 2, x[i, ]))
    both <- sweep(x, 2, x[i, ], "+")
    d[i, ] <- rowSums(diff) / rowSums(both)
  }
  as.dist(d)
}

dissimilarity <- bray_curtis(profile_matrix)

pcoa <- cmdscale(dissimilarity, k = 2, eig = TRUE)
explained <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]))

ordination <- tibble(sample = rownames(pcoa$points),
                     PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2]) %>%
  inner_join(samples, by = c("sample" = "sample_metag"))

ggplot(ordination, aes(PC1, PC2, colour = Layer)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = paste0("PCo1 (", explained[1], "%)"), y = paste0("PCo2 (", explained[2], "%)"))
```

**[Step 5: a biological question](#step-5-a-biological-question)**

```r
cspa <- read_csv("data/tara_cspa.csv", show_col_types = FALSE)

cspa_genes <- samples %>%
  inner_join(cspa, by = c("sample_metag" = "sample")) %>%
  rename(genes_per_cell = cspa_per_cell)

ggplot(cspa_genes, aes(Temperature, genes_per_cell, colour = polar)) +
  geom_point(alpha = 0.8) +
  scale_y_log10() +
  labs(x = "Temperature (°C)", y = "cspA genes per cell (log scale)")

cor.test(cspa_genes$genes_per_cell, cspa_genes$Temperature, method = "spearman")
```

:::
