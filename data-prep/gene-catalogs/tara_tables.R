# Convert the Tara Oceans gene-catalog subset (Salazar et al. 2019, Cell) into compact
# long-format tables that WebR can load quickly in the browser.
# Usage: Rscript tara_tables.R <metat_tutorial_dir> <out_dir>
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]
out_dir <- args[2]

to_long <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    pivot_longer(-c(reference, length, KO, Description),
                 names_to = "sample", values_to = "value") %>%
    filter(value > 0)
}

# Raw insert counts for the 10 marker genes and rpoA/rpoB (K03040/K03043)
mgs <- to_long(file.path(in_dir, "MGS_K03040_K03043_tara.tsv.gz"))

kos <- mgs %>% distinct(KO, Description) %>% arrange(KO)

# 1. Gene-level slice for the hands-on normalisation: 12 metagenomes, all genes
set.seed(1)
samples <- read_csv(file.path(in_dir, "sample_info.csv"), show_col_types = FALSE)
slice_samples <- samples %>%
  group_by(Layer) %>% slice_sample(n = 4) %>% ungroup() %>% pull(sample_metag)
slice <- mgs %>%
  filter(sample %in% slice_samples) %>%
  transmute(gene = reference, length, KO, sample, inserts = as.integer(round(value)))

# 2. Per-sample, per-KO sums for all 258 samples: raw inserts and length-normalised inserts
ko_sums <- mgs %>%
  group_by(sample, KO) %>%
  summarise(inserts = sum(value), inserts_per_bp = sum(value / length), .groups = "drop")

# 3. The 2,000 marker genes found in the most metagenomes, for ordination (metaG only)
metag <- mgs %>% filter(grepl("_METAG$", sample))
top_genes <- metag %>%
  count(reference, name = "n_samples") %>%
  slice_max(n_samples, n = 2000, with_ties = FALSE) %>%
  pull(reference)
top <- metag %>%
  filter(reference %in% top_genes) %>%
  transmute(gene = reference, length, KO, sample, inserts = as.integer(round(value)))

# cspA (K03704), already normalised by gene length and per cell
cspa <- to_long(file.path(in_dir, "K03704_tara_lengthnorm_percell.tsv.gz")) %>%
  group_by(sample) %>%
  summarise(cspa_per_cell = sum(value), .groups = "drop")

write_csv(slice, file.path(out_dir, "tara_gene_slice.csv.gz"))
write_csv(ko_sums, file.path(out_dir, "tara_ko_sums.csv"))
write_csv(top, file.path(out_dir, "tara_top_genes.csv.gz"))
write_csv(kos, file.path(out_dir, "tara_kos.csv"))
write_csv(cspa, file.path(out_dir, "tara_cspa.csv"))
write_csv(samples, file.path(out_dir, "tara_samples.csv"))

cat(sprintf("slice: %d rows; KO sums: %d rows; top genes: %d rows\n",
            nrow(slice), nrow(ko_sums), nrow(top)))
