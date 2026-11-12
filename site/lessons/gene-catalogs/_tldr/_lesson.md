## The whole pipeline at a glance

The commands of every part, in order, without the explanations. Each part has its own summary at its top, with links to the sections.

::: {.callout-tip collapse="true" title="Show all the commands"}

**Terminal (bash)**

```bash
# ---- A. Gene calling: How were the assemblies made?
grep ">" data/scaffolds/metag1.scaffolds.min500.fasta | head -3
grep -c ">" data/scaffolds/*.fasta

# ---- A. Gene calling: Predict genes with Prodigal
mkdir -p genes
for s in metag1 metag2 metag3; do
  prodigal -p meta -c -q \
    -i data/scaffolds/$s.scaffolds.min500.fasta \
    -a genes/$s.faa -d genes/$s.fna -f gff -o genes/$s.gff
done
ls genes

# ---- A. Gene calling: Look at the output
grep -c ">" genes/*.faa

# ---- B. Clustering: Build the catalog with MMseqs2
mkdir -p catalog
cat genes/*.fna | gzip > catalog/genes.fna.gz
cat genes/*.faa | gzip > catalog/genes.faa.gz
cd catalog

mmseqs createdb genes.fna.gz genes.db --dbtype 2 --shuffle 0

mmseqs cluster genes.db genes.cluster tmp \
    --kmer-per-seq-scale 0 \
    --kmer-per-seq 1000 \
    -s 4 \
    --max-seq-len 80000 \
    --remove-tmp-files 0 \
    --cluster-mode 2 \
    --min-seq-id 0.95 \
    --threads 4 \
    --cov-mode 1 \
    -c 0.9 \
    --spaced-kmer-mode 0 \
    --alignment-mode 3 \
    --cluster-reassign 1

mmseqs createtsv genes.db genes.db genes.cluster genes.cluster.tsv
mmseqs version > mmseqs.version
cd ..

# ---- B. Clustering: Name the clusters and write the catalog
Rscript solutions/name_clusters.R catalog/genes.cluster.tsv catalog/genes.fna.gz catalog/clusters.tsv
cut -f 2 catalog/clusters.tsv | tail -n +2 > catalog/representatives.txt
seqtk subseq catalog/genes.fna.gz catalog/representatives.txt > catalog/representatives.fna
seqtk subseq catalog/genes.faa.gz catalog/representatives.txt > catalog/representatives.faa
for ext in fna faa; do
  awk 'NR==FNR { if (FNR > 1) name[$2] = $1; next }
       /^>/ { id = substr($1, 2); print ">" name[id]; next } { print }' \
      catalog/clusters.tsv catalog/representatives.$ext > catalog/catalog.$ext
done
head -4 catalog/clusters.tsv

# ---- C. Quantification: Map the reads
mkdir -p alignments
bwa index catalog/catalog.fna

for s in metag1 metag2 metag3; do
  for r in 1 2; do
    bwa mem -a -t 4 catalog/catalog.fna data/reads/$s.$r.fq.gz \
      | samtools view -F 4 -b - > alignments/$s.r$r.bam
  done
done

# ---- C. Quantification: Pull out what we need
for s in metag1 metag2 metag3; do
  for r in 1 2; do
    samtools view alignments/$s.r$r.bam \
      | awk -v OFS='\t' '{ nm = "NA"; for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6);
                           print $1, $2, $3, $6, nm }' \
      > alignments/$s.r$r.tsv
  done
done
head -6 alignments/metag1.r1.tsv

# ---- C. Quantification: Run it for all samples
Rscript solutions/count_genes.R alignments data/reads counts
head counts/gene_counts.tsv
```

**R**

```r
# ---- C. Quantification: Filter: keep only confident alignments
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())
alignments <- read_csv("data/alignments_metag1_sample.csv", show_col_types = FALSE)

aligned_length <- function(cigar) {
  ops <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
  n <- as.integer(sub("[A-Z=]$", "", ops))
  op <- sub("^[0-9]+", "", ops)
  sum(n[op %in% c("M", "I", "D", "=", "X")])
}

alignments <- alignments %>%
  mutate(aln_length = vapply(cigar, aligned_length, numeric(1)),
         identity = (aln_length - nm) / aln_length,
         passed = aln_length >= 45 & identity >= 0.95)

alignments %>% count(passed) %>% mutate(fraction = n / sum(n))

ggplot(alignments, aes(aln_length, identity, colour = passed)) +
  geom_jitter(alpha = 0.4, width = 0.5, height = 0.002) +
  geom_vline(xintercept = 45, linetype = 2) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  labs(x = "Aligned length (bp)", y = "Identity")

# ---- C. Quantification: Count inserts per gene
inserts <- alignments %>%
  filter(passed) %>%
  group_by(insert) %>%
  summarise(n_genes = n_distinct(gene), gene = first(gene), .groups = "drop")

inserts %>% count(n_genes)

# ---- C. Quantification: What about multi-mapping inserts?
inserts %>%
  filter(n_genes == 1) %>%
  count(gene, sort = TRUE)

# ---- D. Downstream analysis: The data
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())

samples <- read_csv("data/tara_samples.csv", show_col_types = FALSE)
kos <- read_csv("data/tara_kos.csv", show_col_types = FALSE)
genes <- read_csv("data/tara_gene_slice.csv.gz", show_col_types = FALSE)
cat("Loaded", nrow(samples), "samples and", nrow(genes), "gene counts.\n")

# ---- D. Downstream analysis: Why divide before summing?
ko_sums <- read_csv("data/tara_ko_sums.csv", show_col_types = FALSE)

ko_sums %>%
  filter(KO %in% c("K03040", "K03043"), grepl("_METAG$", sample)) %>%
  select(sample, KO, inserts, inserts_per_bp) %>%
  pivot_wider(names_from = KO, values_from = c(inserts, inserts_per_bp)) %>%
  summarise(
    raw_ratio = median(inserts_K03043 / inserts_K03040),
    length_normalised_ratio = median(inserts_per_bp_K03043 / inserts_per_bp_K03040)
  )

# ---- D. Downstream analysis: Step 2: sequencing depth and cells
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

# ---- D. Downstream analysis: Why the median of 10 genes?
per_cell <- ko_sums %>%
  inner_join(mg_median, by = "sample") %>%
  mutate(per_cell = inserts_per_bp / mg_median)

per_cell %>%
  filter(grepl("_METAG$", sample)) %>%
  group_by(KO) %>%
  summarise(median_copies_per_cell = median(per_cell)) %>%
  left_join(kos, by = "KO")

# ---- D. Downstream analysis: Step 4: comparing samples
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

# ---- D. Downstream analysis: Step 5: a biological question
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

Parts: [A. Gene calling](01-gene-calling.qmd) · [B. Clustering](02-clustering.qmd) · [C. Quantification](03-quantification.qmd) · [D. Downstream analysis](04-downstream.qmd)
