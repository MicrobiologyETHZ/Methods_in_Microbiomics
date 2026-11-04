# Build the data files for the metatranscriptomics lesson (site/lessons/metatranscriptomics/data).
#
#   Rscript data-prep/metatranscriptomics/prepare.R
#
# Input: docs/downloads/metat_defined_community_data.tar.gz (OMM12 defined community, LPS vs PBS).
# KEGG names are looked up once with KEGGREST and cached in data-prep/metatranscriptomics/kegg_*.csv,
# so later runs work offline. Delete the cache files to refresh them.
# The Tara files for Part 5 come from the gene-catalog lesson (see the end of this script).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(DESeq2)
})

here <- "data-prep/metatranscriptomics"
out <- "site/lessons/metatranscriptomics/data"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

tmp <- tempfile()
untar("docs/downloads/metat_defined_community_data.tar.gz", exdir = tmp, tar = "internal")
load(file.path(tmp, "metat_defined_community_data", "METAT_raw_data.RData"))

na_dash <- function(x) if_else(x %in% c("-", "nan", ""), NA_character_, x)

# --- Samples, genes, counts -------------------------------------------------------------------
samples <- sample_data %>%
  select(sample_id, Treatment) %>%
  arrange(Treatment, sample_id)
write_csv(samples, file.path(out, "samples.csv"))

genes <- annotations %>%
  transmute(
    locus_tag, genome,
    length = End - Start + 1,
    KO = na_dash(KO),
    # keep the KEGG "ko" pathway IDs only; the "map" IDs are duplicates of them
    pathways = na_dash(KEGG_Pathway) %>%
      strsplit(",") %>%
      vapply(function(p) paste(p[startsWith(p, "ko")], collapse = ";"), "") %>%
      na_if("")
  )

counts <- raw_counts %>%
  select(locus_tag, all_of(samples$sample_id)) %>%
  mutate(across(-locus_tag, as.integer))
stopifnot(identical(counts$locus_tag, genes$locus_tag))
write_csv(counts, file.path(out, "counts.csv.gz"))

# --- KEGG names (cached) ----------------------------------------------------------------------
cache <- function(file, fetch) {
  path <- file.path(here, file)
  if (!file.exists(path)) write_csv(fetch(), path)
  read_csv(path, show_col_types = FALSE)
}

ko_names <- cache("kegg_ko.csv", function() {
  x <- KEGGREST::keggList("ko")
  tibble(KO = sub("^ko:", "", names(x)), name = unname(x))
})
pathway_names <- cache("kegg_pathways.csv", function() {
  x <- KEGGREST::keggList("pathway", "ko")
  tibble(pathway = sub("^path:", "", names(x)), name = unname(x))
})

used_kos <- unique(unlist(strsplit(na.omit(genes$KO), ",")))
write_csv(filter(ko_names, KO %in% used_kos), file.path(out, "kegg_ko.csv"))
# drop pathway IDs that KEGG has retired (they have no name any more)
genes <- genes %>%
  mutate(pathways = strsplit(coalesce(pathways, ""), ";") %>%
           vapply(function(p) paste(p[p %in% pathway_names$pathway], collapse = ";"), "") %>%
           na_if(""))
write_csv(genes, file.path(out, "genes.csv.gz"))
used_pathways <- unique(unlist(strsplit(na.omit(genes$pathways), ";")))
write_csv(filter(pathway_names, pathway %in% used_pathways), file.path(out, "kegg_pathways.csv"))

# --- Differential expression results for Part 4 ---------------------------------------------------
# The method taught in Part 3: one DESeq2 model, with size factors estimated within each genome
# and passed in as normalisation factors (taxon-specific scaling). Part 4 uses I48 only.
de_genomes <- c("I48", "YL32", "YL58")

m <- counts %>%
  filter(genes$genome %in% de_genomes) %>%
  tibble::column_to_rownames("locus_tag") %>%
  as.matrix()
genome <- genes$genome[match(rownames(m), genes$locus_tag)]

nf <- m
for (g in de_genomes) {
  rows <- genome == g
  sf <- estimateSizeFactorsForMatrix(m[rows, ])
  nf[rows, ] <- matrix(sf, sum(rows), ncol(m), byrow = TRUE)
}

coldata <- samples %>%
  mutate(Treatment = factor(Treatment, levels = c("PBS", "LPS"))) %>%
  tibble::column_to_rownames("sample_id")
dds <- DESeqDataSetFromMatrix(m, coldata[colnames(m), , drop = FALSE], design = ~ Treatment)
normalizationFactors(dds) <- nf
dds <- DESeq(dds, quiet = TRUE)

de <- results(dds, contrast = c("Treatment", "LPS", "PBS"), alpha = 0.01) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("locus_tag") %>%
  mutate(genome = genome[match(locus_tag, rownames(m))], .after = locus_tag) %>%
  select(locus_tag, genome, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
write_csv(filter(de, genome == "I48"), file.path(out, "de_i48.csv.gz"))

# --- Tara (Part 5) ------------------------------------------------------------------------------
# Built by data-prep/gene-catalogs/tara_tables.R; copied so each lesson's data folder is complete.
gc <- "site/lessons/gene-catalogs/data"
for (f in c("tara_samples.csv", "tara_cspa.csv", "tara_ko_sums.csv", "tara_kos.csv")) {
  file.copy(file.path(gc, f), file.path(out, f), overwrite = TRUE)
}

jsonlite::write_json(
  list(
    source = "docs/downloads/metat_defined_community_data.tar.gz",
    deseq2_version = as.character(packageVersion("DESeq2")),
    de_genomes = de_genomes
  ),
  file.path(out, "provenance.json"), auto_unbox = TRUE, pretty = TRUE
)

cat("Wrote", out, "\n")
print(file.info(list.files(out, full.names = TRUE))[, "size", drop = FALSE])
