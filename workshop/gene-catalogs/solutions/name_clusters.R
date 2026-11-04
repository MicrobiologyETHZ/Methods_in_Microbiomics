# Summarise the MMseqs2 cluster table and name clusters the way mCATs does:
# catalog_000000001 is the cluster with the longest representative.
# Usage: Rscript name_clusters.R genes.cluster.tsv genes.fna.gz clusters.tsv
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
cluster_file <- args[1]
fasta_file <- args[2]
out_file <- args[3]

# gene lengths from the FASTA file (sequence length per header)
read_lengths <- function(path) {
  lines <- readLines(path)
  is_header <- startsWith(lines, ">")
  ids <- sub(" .*", "", sub("^>", "", lines[is_header]))
  seq_number <- cumsum(is_header)
  lengths <- tapply(nchar(lines[!is_header]), seq_number[!is_header], sum)
  tibble(gene = ids, length = as.integer(lengths))
}

gene_lengths <- read_lengths(fasta_file)

clusters <- read.delim(cluster_file, header = FALSE,
                       col.names = c("representative", "member"))

summary <- clusters %>%
  count(representative, name = "size") %>%
  left_join(gene_lengths, by = c("representative" = "gene")) %>%
  rename(rep_length = length) %>%
  arrange(desc(rep_length), representative) %>%
  mutate(cluster_id = sprintf("catalog_%09d", row_number())) %>%
  select(cluster_id, representative, rep_length, size)

write.table(summary, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d genes in %d clusters\n", nrow(clusters), nrow(summary)))
