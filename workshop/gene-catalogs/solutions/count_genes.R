# Filter alignments and count uniquely mapped inserts per gene.
# Usage: Rscript count_genes.R <alignments_dir> <reads_dir> <out_dir>
suppressMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
aln_dir <- args[1]
reads_dir <- args[2]
out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE)

MIN_ALN_LENGTH <- 45
MIN_IDENTITY <- 0.95

# Number of alignment columns in a CIGAR string: matches/mismatches (M, =, X),
# insertions (I) and deletions (D). Clipped bases (S, H) are not aligned.
cigar_aln_length <- function(cigar) {
  ops <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))
  vapply(ops, function(o) {
    n <- as.integer(sub("[A-Z=]$", "", o))
    op <- sub("^[0-9]+", "", o)
    sum(n[op %in% c("M", "I", "D", "=", "X")])
  }, numeric(1))
}

read_alignments <- function(file) {
  read.delim(file, header = FALSE,
             col.names = c("insert", "flag", "gene", "cigar", "nm"),
             colClasses = c("character", "integer", "character", "character", "integer"))
}

samples <- sub("\\.r1\\.tsv$", "", list.files(aln_dir, pattern = "\\.r1\\.tsv$"))
counts <- list()
stats <- list()

for (s in samples) {
  aln <- bind_rows(
    read_alignments(file.path(aln_dir, paste0(s, ".r1.tsv"))),
    read_alignments(file.path(aln_dir, paste0(s, ".r2.tsv")))
  ) %>%
    mutate(aln_length = cigar_aln_length(cigar),
           identity = (aln_length - nm) / aln_length)

  passed <- aln %>% filter(aln_length >= MIN_ALN_LENGTH, identity >= MIN_IDENTITY)

  # An insert (read pair) is unique if all its passing alignments, from either mate, hit one gene
  inserts <- passed %>%
    group_by(insert) %>%
    summarise(n_genes = n_distinct(gene), gene = first(gene), .groups = "drop")

  counts[[s]] <- inserts %>%
    filter(n_genes == 1) %>%
    count(gene, name = s)

  total_inserts <- length(readLines(file.path(reads_dir, paste0(s, ".1.fq.gz")))) / 4
  stats[[s]] <- tibble(
    sample = s,
    total_inserts = total_inserts,
    mapped_inserts = nrow(inserts),
    unique_inserts = sum(inserts$n_genes == 1),
    multi_inserts = sum(inserts$n_genes > 1),
    mapping_rate = round(nrow(inserts) / total_inserts, 4)
  )
}

count_table <- Reduce(function(a, b) full_join(a, b, by = "gene"), counts) %>%
  mutate(across(-gene, ~ replace_na(.x, 0L))) %>%
  arrange(gene)

write.table(count_table, file.path(out_dir, "gene_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(stats), file.path(out_dir, "mapping_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(bind_rows(stats))
