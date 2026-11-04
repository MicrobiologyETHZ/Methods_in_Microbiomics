# Summarise a finished lesson run into the metrics that `check` compares against.
# Usage: Rscript reference_summary.R <run_dir> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "solutions", "metrics.R"))
write.table(lesson_metrics(args[1]), args[2], sep = "\t", quote = FALSE, row.names = FALSE)
