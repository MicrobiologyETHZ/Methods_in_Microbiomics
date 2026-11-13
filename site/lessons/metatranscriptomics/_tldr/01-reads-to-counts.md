::: {.callout-note .callout-glance appearance="simple" collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**[The workflow](#the-workflow)**

```bash
cat genomes/*.fna > community.fna
bowtie2-build --threads 4 community.fna community

bowtie2 -p 4 -x community -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  | samtools sort -@ 4 -o sample.bam

cat annotations/*.gff > community.gff

featureCounts -T 4 -p --countReadPairs -s 2 -Q 10 \
  -t CDS -g locus_tag \
  -a community.gff -o counts.tsv sample*.bam
```

:::
