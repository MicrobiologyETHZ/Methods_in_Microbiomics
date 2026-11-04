# Metrics describing a lesson run, shared by `check` and the data-prep reference summary.
# Each metric belongs to a part (a, b or c) and is NA when its file does not exist yet.

count_fasta <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  on.exit(close(con))
  sum(startsWith(readLines(con), ">"))
}

lesson_metrics <- function(run_dir) {
  p <- function(...) file.path(run_dir, ...)
  samples <- c("metag1", "metag2", "metag3")
  rows <- list()
  add <- function(part, metric, value) {
    rows[[length(rows) + 1]] <<- data.frame(part = part, metric = metric, value = value)
  }

  for (s in samples) add("a", paste0(s, "_genes"), count_fasta(p("genes", paste0(s, ".faa"))))

  cluster_tsv <- p("catalog", "genes.cluster.tsv")
  if (file.exists(cluster_tsv)) {
    clusters <- read.delim(cluster_tsv, header = FALSE)
    add("b", "genes_clustered", nrow(clusters))
    add("b", "clusters", length(unique(clusters$V1)))
  } else {
    add("b", "genes_clustered", NA)
    add("b", "clusters", NA)
  }
  add("b", "catalog_sequences", count_fasta(p("catalog", "catalog.fna")))

  counts_file <- p("counts", "gene_counts.tsv")
  if (file.exists(counts_file)) {
    counts <- read.delim(counts_file)
    for (s in samples) add("c", paste0(s, "_unique_inserts"), sum(counts[[s]]))
    add("c", "genes_detected", nrow(counts))
  } else {
    for (s in samples) add("c", paste0(s, "_unique_inserts"), NA)
    add("c", "genes_detected", NA)
  }
  do.call(rbind, rows)
}
