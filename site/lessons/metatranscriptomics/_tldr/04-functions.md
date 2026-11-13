::: {.callout-note .callout-glance appearance="simple" collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**[The data](#the-data)**

```r
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
cat("Ready: packages and data are loaded.\n")

universe <- unique(pathway_genes$locus_tag)
down <- i48 %>% filter(padj < 0.01, log2FoldChange < -1) %>% pull(locus_tag) %>% intersect(universe)
up <- i48 %>% filter(padj < 0.01, log2FoldChange > 1) %>% pull(locus_tag) %>% intersect(universe)
cat(length(universe), "genes in the universe;", length(down), "down and", length(up), "up in LPS mice\n")
```

**[Over-representation analysis](#over-representation-analysis)**

```r
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
```

**[Gene set enrichment analysis](#gene-set-enrichment-analysis)**

```r
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
```

:::
