## The whole pipeline at a glance

The commands of every part, in order, without the explanations. Each part has its own summary at its top, with links to the sections.

::: {.callout-tip collapse="true" title="Show all the commands"}

**R**

```r
# ---- 2. Explore the data: The data
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

# ---- 2. Explore the data: Step 1: who dominates the transcriptome?
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

# ---- 2. Explore the data: Step 2: is there enough data for each strain?
depth <- strain_inserts %>%
  group_by(genome, Treatment) %>%
  summarise(n_adequate = sum(inserts >= 5e5), .groups = "drop")

ggplot(depth, aes(Treatment, genome, fill = n_adequate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_adequate)) +
  scale_fill_gradient(low = "grey95", high = "steelblue", limits = c(0, 6)) +
  labs(x = NULL, y = NULL, fill = "Samples with\n≥ 500k inserts")

# ---- 2. Explore the data: Step 3: normalising for exploration
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

# ---- 2. Explore the data: Step 4: do replicates agree?
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

# ---- 2. Explore the data: Step 5: how fast are the strains growing?
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

# ---- 3. Differential expression: The data
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

# ---- 3. Differential expression: Step 1: size factors
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

# ---- 3. Differential expression: Step 2: DESeq2 with per-strain size factors
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

# ---- 3. Differential expression: Step 3: look at the results
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

# ---- 4. From genes to functions: The data
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())

genes <- read_csv("data/genes.csv.gz", show_col_types = FALSE)
de <- read_csv("data/de_i48.csv.gz", show_col_types = FALSE)
pathway_names <- read_csv("data/kegg_pathways.csv", show_col_types = FALSE)

# genes that were tested (padj is NA for genes DESeq2 filtered out, see Part 3)
i48 <- de %>%
  filter(!is.na(padj)) %>%
  left_join(select(genes, locus_tag, pathways), by = "locus_tag")

# one row per gene and pathway
pathway_genes <- i48 %>%
  filter(!is.na(pathways)) %>%
  separate_rows(pathways, sep = ";") %>%
  distinct(pathway = pathways, locus_tag)

cat(nrow(i48), "I48 genes tested;", n_distinct(pathway_genes$locus_tag),
    "of them are in at least one KEGG pathway\n")

universe <- unique(pathway_genes$locus_tag)
down <- i48 %>% filter(padj < 0.01, log2FoldChange < -1) %>% pull(locus_tag) %>% intersect(universe)
up <- i48 %>% filter(padj < 0.01, log2FoldChange > 1) %>% pull(locus_tag) %>% intersect(universe)
cat(length(universe), "genes in the universe;", length(down), "down and", length(up), "up in LPS mice\n")

# ---- 4. From genes to functions: Over-representation analysis
ora <- function(selected) {
  pathway_genes %>%
    group_by(pathway) %>%
    summarise(size = n(), hits = sum(locus_tag %in% selected)) %>%
    filter(size >= 10) %>%
    mutate(
      expected = size * length(selected) / length(universe),
      ratio = hits / expected,
      p = phyper(hits - 1, size, length(universe) - size, length(selected), lower.tail = FALSE),
      padj = p.adjust(p, method = "BH")
    ) %>%
    left_join(pathway_names, by = "pathway") %>%
    arrange(p) %>%
    select(pathway, name, size, hits, expected, ratio, padj)
}

ora(down) %>% head(8)
ora(up) %>% head(8)

# ---- 4. From genes to functions: Gene set enrichment analysis
library(fgsea)

ranks <- i48 %>%
  filter(locus_tag %in% universe) %>%
  arrange(desc(stat)) %>%
  { setNames(.$stat, .$locus_tag) }

set.seed(42)
gsea <- fgsea(split(pathway_genes$locus_tag, pathway_genes$pathway), ranks,
              minSize = 10, maxSize = 500,
              BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))

gsea %>%
  as_tibble() %>%
  left_join(pathway_names, by = "pathway") %>%
  filter(padj < 0.05) %>%
  arrange(NES) %>%
  select(pathway, name, size, NES, padj)

# ---- 5. When you have metagenomes: From strains to whole communities
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())

samples <- read_csv("data/tara_samples.csv", show_col_types = FALSE)
kos <- read_csv("data/tara_kos.csv", show_col_types = FALSE) %>%
  add_row(KO = "K03704", Description = "cspA; cold shock protein")
marker_genes <- c("K06942", "K01889", "K01887", "K01875", "K01883",
                  "K01869", "K01873", "K01409", "K03106", "K03110")

# copies (metagenomes) or transcripts (metatranscriptomes) per cell, for 12 KOs, plus cspA
per_cell <- read_csv("data/tara_ko_sums.csv", show_col_types = FALSE) %>%
  group_by(sample) %>%
  mutate(per_cell = inserts_per_bp / median(inserts_per_bp[KO %in% marker_genes])) %>%
  ungroup() %>%
  select(sample, KO, per_cell) %>%
  bind_rows(read_csv("data/tara_cspa.csv", show_col_types = FALSE) %>%
              transmute(sample, KO = "K03704", per_cell = cspa_per_cell))

cat(n_distinct(per_cell$sample), "samples,", n_distinct(per_cell$KO), "KOs\n")

# ---- 5. When you have metagenomes: Cold-shock protein: genes or expression?
cspa_genes <- samples %>%
  inner_join(filter(per_cell, KO == "K03704"), by = c("sample_metag" = "sample")) %>%
  rename(genes_per_cell = per_cell)

ggplot(cspa_genes, aes(Temperature, genes_per_cell, colour = polar)) +
  geom_point(alpha = 0.8) +
  scale_y_log10() +
  labs(x = "Temperature (°C)", y = "cspA genes per cell (log scale)")

cor.test(cspa_genes$genes_per_cell, cspa_genes$Temperature, method = "spearman", exact = FALSE)

cspa_both <- cspa_genes %>%
  inner_join(filter(per_cell, KO == "K03704"), by = c("sample_metat" = "sample", "KO")) %>%
  rename(transcripts_per_cell = per_cell) %>%
  mutate(expression = transcripts_per_cell / genes_per_cell)

ggplot(cspa_both, aes(Temperature, expression, colour = polar)) +
  geom_point(alpha = 0.8) +
  scale_y_log10() +
  labs(x = "Temperature (°C)", y = "cspA expression per gene copy (log scale)")

cor.test(cspa_both$expression, cspa_both$Temperature, method = "spearman", exact = FALSE)

# ---- 5. When you have metagenomes: How much is turnover, how much is regulation?
paired <- samples %>%
  select(sample_metag, sample_metat, Temperature) %>%
  inner_join(per_cell, by = c("sample_metag" = "sample")) %>%
  rename(genes = per_cell) %>%
  inner_join(per_cell, by = c("sample_metat" = "sample", "KO")) %>%
  rename(transcripts = per_cell) %>%
  filter(genes > 0, transcripts > 0) %>%
  mutate(expression = transcripts / genes)

variance <- paired %>%
  group_by(KO) %>%
  summarise(
    # covariance first: below, `genes` and `expression` are replaced by their variances
    covariance = 2 * cov(log2(genes), log2(expression)),
    genes = var(log2(genes)),
    expression = var(log2(expression))
  ) %>%
  left_join(kos, by = "KO") %>%
  mutate(label = paste(KO, sub(";.*", "", Description)))

variance %>%
  pivot_longer(c(genes, expression, covariance), names_to = "source", values_to = "variance") %>%
  mutate(source = factor(source, levels = c("covariance", "expression", "genes"))) %>%
  ggplot(aes(variance, reorder(label, variance), fill = source)) +
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(genes = "#2f6f8f", expression = "#e8590c", covariance = "grey70")) +
  labs(x = "Variance of log2 transcripts per cell across samples", y = NULL, fill = "Part due to")

shares <- variance %>%
  mutate(total = genes + expression + covariance,
         expression_share = expression / total) %>%
  select(label, total, expression_share) %>%
  arrange(desc(total))
shares
```

:::

Parts: [2. Explore the data](02-explore.qmd) · [3. Differential expression](03-differential-expression.qmd) · [4. From genes to functions](04-functions.qmd) · [5. When you have metagenomes](05-with-metagenomes.qmd)
