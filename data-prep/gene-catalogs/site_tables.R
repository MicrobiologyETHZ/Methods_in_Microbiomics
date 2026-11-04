# Small tables and text excerpts from the reference run, used by the lesson pages.
# Usage: Rscript site_tables.R <run_dir> <out_dir>
suppressMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
run <- args[1]
out <- args[2]
excerpts <- file.path(out, "_excerpts")  # leading _ : Quarto does not render these
dir.create(excerpts, showWarnings = FALSE, recursive = TRUE)
samples <- c("metag1", "metag2", "metag3")

# Write the first lines of a file as a fenced block that pages can include
excerpt <- function(path, name, n = 6, width = 110) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  lines <- readLines(con, n = n)
  close(con)
  lines <- ifelse(nchar(lines) > width, paste0(substr(lines, 1, width), "…"), lines)
  writeLines(c("```text", lines, "```"), file.path(excerpts, paste0(name, ".md")))
}

# Part A: prodigal headers carry coordinates, strand and partial/start information
genes <- bind_rows(lapply(samples, function(s) {
  h <- grep("^>", readLines(file.path(run, "genes", paste0(s, ".faa"))), value = TRUE)
  f <- strsplit(sub("^>", "", h), " # ")
  tibble(
    sample = s,
    gene = vapply(f, `[`, "", 1),
    start = as.integer(vapply(f, `[`, "", 2)),
    end = as.integer(vapply(f, `[`, "", 3)),
    partial = sub(".*partial=([0-9]+).*", "\\1", h),
    start_type = sub(".*start_type=([A-Za-z]+).*", "\\1", h)
  ) %>% mutate(length_nt = end - start + 1)
}))
write_csv(genes, file.path(out, "genes.csv"))
excerpt(file.path(run, "genes", "metag1.gff"), "gff", n = 5)
excerpt(file.path(run, "genes", "metag1.faa"), "faa", n = 4)

# Part B
invisible(file.copy(file.path(run, "catalog", "clusters.tsv"), file.path(out, "clusters.tsv"), overwrite = TRUE))
invisible(file.copy(file.path(run, "catalog", "genes.cluster.tsv"), file.path(out, "genes.cluster.tsv"), overwrite = TRUE))
excerpt(file.path(run, "catalog", "genes.cluster.tsv"), "cluster_tsv", n = 8)
excerpt(file.path(run, "catalog", "clusters.tsv"), "clusters", n = 6)
invisible(file.copy(file.path(run, "catalog", "mmseqs.version"), file.path(out, "mmseqs.version"), overwrite = TRUE))

# Part C: all raw alignments (before filtering) of 3,000 random inserts from metag1, both mates
read_aln <- function(mate) {
  read.delim(file.path(run, "alignments", paste0("metag1.r", mate, ".tsv")), header = FALSE,
             col.names = c("insert", "flag", "gene", "cigar", "nm")) %>%
    mutate(mate = mate, .after = insert)
}
aln <- bind_rows(read_aln(1), read_aln(2))
set.seed(1)
picked <- sample(unique(aln$insert), 3000)
write_csv(aln %>% filter(insert %in% picked) %>% arrange(insert, mate),
          file.path(out, "alignments_metag1_sample.csv"))
excerpt(file.path(run, "alignments", "metag1.r1.tsv"), "alignments", n = 6)
invisible(file.copy(file.path(run, "counts", "gene_counts.tsv"), file.path(out, "gene_counts.tsv"), overwrite = TRUE))
invisible(file.copy(file.path(run, "counts", "mapping_stats.tsv"), file.path(out, "mapping_stats.tsv"), overwrite = TRUE))
