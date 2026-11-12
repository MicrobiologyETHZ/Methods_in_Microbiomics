::: {.callout-tip collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**The data** [↗](#the-data)

```bash
metaspades.py -t 4 -m 12 --only-assembler \
  --pe1-1 reads/metag1.1.fq.gz --pe1-2 reads/metag1.2.fq.gz -o asm/metag1
python scaffold_filter.py metag1 scaffolds asm/metag1/scaffolds.fasta asm/metag1 META

grep ">" data/scaffolds/metag1.scaffolds.min500.fasta | head -3
grep -c ">" data/scaffolds/*.fasta
```

**Predict genes with Prodigal** [↗](#predict-genes-with-prodigal)

```bash
mkdir -p genes
for s in metag1 metag2 metag3; do
  prodigal -p meta -c -q \
    -i data/scaffolds/$s.scaffolds.min500.fasta \
    -a genes/$s.faa -d genes/$s.fna -f gff -o genes/$s.gff
done
ls genes
```

**Look at the output** [↗](#look-at-the-output)

```bash
grep -c ">" genes/*.faa
```

:::
