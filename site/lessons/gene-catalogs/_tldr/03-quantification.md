::: {.callout-tip collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**Map the reads** [↗](#map-the-reads)

```bash
mkdir -p alignments
bwa index catalog/catalog.fna

for s in metag1 metag2 metag3; do
  for r in 1 2; do
    bwa mem -a -t 4 catalog/catalog.fna data/reads/$s.$r.fq.gz \
      | samtools view -F 4 -b - > alignments/$s.r$r.bam
  done
done
```

**Pull out what we need** [↗](#pull-out-what-we-need)

```bash
for s in metag1 metag2 metag3; do
  for r in 1 2; do
    samtools view alignments/$s.r$r.bam \
      | awk -v OFS='\t' '{ nm = "NA"; for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6);
                           print $1, $2, $3, $6, nm }' \
      > alignments/$s.r$r.tsv
  done
done
head -6 alignments/metag1.r1.tsv
```

**Filter: keep only confident alignments** [↗](#filter-keep-only-confident-alignments)

```r
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
```

**Count inserts per gene** [↗](#count-inserts-per-gene)

```r
inserts <- alignments %>%
  filter(passed) %>%
  group_by(insert) %>%
  summarise(n_genes = n_distinct(gene), gene = first(gene), .groups = "drop")

inserts %>% count(n_genes)
```

**What about multi-mapping inserts?** [↗](#what-about-multi-mapping-inserts)

```r
inserts %>%
  filter(n_genes == 1) %>%
  count(gene, sort = TRUE)
```

**Run it for all samples** [↗](#run-it-for-all-samples)

```bash
Rscript solutions/count_genes.R alignments data/reads counts
head counts/gene_counts.tsv
```

:::
