::: {.callout-note .callout-glance appearance="simple" collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**[From strains to whole communities](#from-strains-to-whole-communities)**

```r
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
theme_set(theme_bw())

samples <- read_csv("data/tara_samples.csv", show_col_types = FALSE)
kos <- read_csv("data/tara_kos.csv", show_col_types = FALSE) %>%
  add_row(KO = "K03704", Description = "cspA; cold shock protein")
marker_genes <- c("K06942", "K01889", "K01887", "K01875", "K01883",
                  "K01869", "K01873", "K01409", "K03106", "K03110")

# copies (metagenomes) or transcripts (metatranscriptomes) per cell, for 12 KOs, plus cspA
per_cell <- read_csv("data/tara_ko_sums.csv", show_col_types = FALSE) %>%
  group_by(sample) %>%
  mutate(per_cell = inserts_per_bp / median(inserts_per_bp[KO %in% marker_genes])) %>%
  ungroup() %>%
  select(sample, KO, per_cell) %>%
  bind_rows(read_csv("data/tara_cspa.csv", show_col_types = FALSE) %>%
              transmute(sample, KO = "K03704", per_cell = cspa_per_cell))

cat(n_distinct(per_cell$sample), "samples,", n_distinct(per_cell$KO), "KOs\n")
```

**[Cold-shock protein: genes or expression?](#cold-shock-protein-genes-or-expression)**

```r
cspa_genes <- samples %>%
  inner_join(filter(per_cell, KO == "K03704"), by = c("sample_metag" = "sample")) %>%
  rename(genes_per_cell = per_cell)

ggplot(cspa_genes, aes(Temperature, genes_per_cell, colour = polar)) +
  geom_point(alpha = 0.8) +
  scale_y_log10() +
  labs(x = "Temperature (°C)", y = "cspA genes per cell (log scale)")

cor.test(cspa_genes$genes_per_cell, cspa_genes$Temperature, method = "spearman", exact = FALSE)

cspa_both <- cspa_genes %>%
  inner_join(filter(per_cell, KO == "K03704"), by = c("sample_metat" = "sample", "KO")) %>%
  rename(transcripts_per_cell = per_cell) %>%
  mutate(expression = transcripts_per_cell / genes_per_cell)

ggplot(cspa_both, aes(Temperature, expression, colour = polar)) +
  geom_point(alpha = 0.8) +
  scale_y_log10() +
  labs(x = "Temperature (°C)", y = "cspA expression per gene copy (log scale)")

cor.test(cspa_both$expression, cspa_both$Temperature, method = "spearman", exact = FALSE)
```

**[How much is turnover, how much is regulation?](#how-much-is-turnover-how-much-is-regulation)**

```r
paired <- samples %>%
  select(sample_metag, sample_metat, Temperature) %>%
  inner_join(per_cell, by = c("sample_metag" = "sample")) %>%
  rename(genes = per_cell) %>%
  inner_join(per_cell, by = c("sample_metat" = "sample", "KO")) %>%
  rename(transcripts = per_cell) %>%
  filter(genes > 0, transcripts > 0) %>%
  mutate(expression = transcripts / genes)

variance <- paired %>%
  group_by(KO) %>%
  summarise(
    # covariance first: below, `genes` and `expression` are replaced by their variances
    covariance = 2 * cov(log2(genes), log2(expression)),
    genes = var(log2(genes)),
    expression = var(log2(expression))
  ) %>%
  left_join(kos, by = "KO") %>%
  mutate(label = paste(KO, sub(";.*", "", Description)))

variance %>%
  pivot_longer(c(genes, expression, covariance), names_to = "source", values_to = "variance") %>%
  mutate(source = factor(source, levels = c("covariance", "expression", "genes"))) %>%
  ggplot(aes(variance, reorder(label, variance), fill = source)) +
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(genes = "#2f6f8f", expression = "#e8590c", covariance = "grey70")) +
  labs(x = "Variance of log2 transcripts per cell across samples", y = NULL, fill = "Part due to")

shares <- variance %>%
  mutate(total = genes + expression + covariance,
         expression_share = expression / total) %>%
  select(label, total, expression_share) %>%
  arrange(desc(total))
shares
```

:::
