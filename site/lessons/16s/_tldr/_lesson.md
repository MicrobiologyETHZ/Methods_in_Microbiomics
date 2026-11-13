## The whole pipeline at a glance

The commands of every part, in order, without the explanations. Each part has its own summary at its top, with links to the sections.

::: {.callout-tip collapse="true" title="Show all the commands"}

**Terminal (bash)**

```bash
# ---- 1. From reads to ASVs: Step 1: primers
FWD=CCTACGGGNGGCWGCAG       # 341F
REV=GACTACHVGGGTATCTAATCC   # 785R (also called 805R)
mkdir -p trimmed
for s in $(cat samples.txt); do
  # keep read pairs that start with the primers, and cut them off
  cutadapt -O 12 --discard-untrimmed --pair-adapters --minimum-length 75 -j 4 \
    -g $FWD -G $REV \
    -o trimmed/$s.tmp_R1.fastq.gz -p trimmed/$s.tmp_R2.fastq.gz \
    raw/${s}_R1.fastq.gz raw/${s}_R2.fastq.gz
  # rarely, a primer was added more than once: remove up to 5 more copies
  cutadapt -O 12 --times 5 -j 4 -g $FWD -o trimmed/$s.tmp2_R1.fastq.gz trimmed/$s.tmp_R1.fastq.gz
  cutadapt -O 12 --times 5 -j 4 -g $REV -o trimmed/$s.tmp2_R2.fastq.gz trimmed/$s.tmp_R2.fastq.gz
  # pair the reads up again, dropping any that became too short
  cutadapt --minimum-length 75 -j 4 \
    -o trimmed/${s}_R1.fastq.gz -p trimmed/${s}_R2.fastq.gz \
    trimmed/$s.tmp2_R1.fastq.gz trimmed/$s.tmp2_R2.fastq.gz
  rm trimmed/$s.tmp*
done
```

**R**

```r
# ---- 1. From reads to ASVs: The reads on this page
library(dada2)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

# forward (R1) and reverse (R2) reads of each sample
path <- system.file("extdata", package = "dada2")
fnFs <- sort(list.files(path, pattern = "^sam.*F\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "^sam.*R\\.fastq\\.gz$", full.names = TRUE))
sample_names <- sub("F\\.fastq\\.gz$", "", basename(fnFs))
sample_names
cat("Ready: packages and data are loaded.\n")

# ---- 1. From reads to ASVs: Step 1: primers
fwd_primer <- "GTGCCAGCMGCCGCGGTAA"   # 515F
rev_primer <- "GGACTACHVGGGTWTCTAAT"  # 806R

# number of reads containing a primer (degenerate bases allowed), in both orientations
primer_hits <- function(primer, fastq) {
  reads <- Biostrings::DNAStringSet(getSequences(fastq))
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer)))
  c(forward = sum(Biostrings::vcountPattern(primer, reads, fixed = FALSE) > 0),
    reverse_complement = sum(Biostrings::vcountPattern(rc, reads, fixed = FALSE) > 0))
}
rbind(fwd_in_R1 = primer_hits(fwd_primer, fnFs[1]),
      rev_in_R1 = primer_hits(rev_primer, fnFs[1]),
      fwd_in_R2 = primer_hits(fwd_primer, fnRs[1]),
      rev_in_R2 = primer_hits(rev_primer, fnRs[1]))

# ---- 1. From reads to ASVs: Step 3: filter and truncate
filtFs <- file.path(tempdir(), "filtered", basename(fnFs))
filtRs <- file.path(tempdir(), "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     rm.phix = TRUE, multithread = FALSE)
out

# ---- 1. From reads to ASVs: Step 4: learn the error rates
errF <- learnErrors(filtFs, multithread = FALSE, verbose = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE, verbose = FALSE)
plotErrors(errF, nominalQ = TRUE)

# ---- 1. From reads to ASVs: Step 5: denoise
ddF <- dada(filtFs, err = errF, multithread = FALSE, verbose = FALSE)
ddR <- dada(filtRs, err = errR, multithread = FALSE, verbose = FALSE)
ddF[[1]]

# ---- 1. From reads to ASVs: Step 6: merge read pairs
mergers <- mergePairs(ddF, filtFs, ddR, filtRs, verbose = FALSE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)                           # samples x ASVs
table(nchar(getSequences(seqtab)))    # ASV lengths

# ---- 1. From reads to ASVs: Step 7: remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = FALSE, verbose = TRUE)
c(ASVs_before = ncol(seqtab), ASVs_after = ncol(seqtab_nochim),
  reads_kept = sum(seqtab_nochim) / sum(seqtab))

# ---- 1. From reads to ASVs: Step 8: how many reads made it?
getN <- function(x) sum(getUniques(x))
track <- data.frame(sample = sample_names,
                    input = out[, "reads.in"], filtered = out[, "reads.out"],
                    denoisedF = sapply(ddF, getN), denoisedR = sapply(ddR, getN),
                    merged = sapply(mergers, getN), nonchim = rowSums(seqtab_nochim),
                    row.names = NULL)
track

track <- track %>%
  mutate(kept = nonchim / input)
track
```

:::

Parts: [1. From reads to ASVs](01-reads-to-asvs.qmd)
