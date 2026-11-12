::: {.callout-tip collapse="true" title="At a glance: the commands of this part"}
For coming back later. Each step is explained in the section it links to.

**Build the catalog with MMseqs2** [↗](#build-the-catalog-with-mmseqs2)

```bash
mkdir -p catalog
cat genes/*.fna | gzip > catalog/genes.fna.gz
cat genes/*.faa | gzip > catalog/genes.faa.gz
cd catalog

mmseqs createdb genes.fna.gz genes.db --dbtype 2 --shuffle 0

mmseqs cluster genes.db genes.cluster tmp \
    --kmer-per-seq-scale 0 \
    --kmer-per-seq 1000 \
    -s 4 \
    --max-seq-len 80000 \
    --remove-tmp-files 0 \
    --cluster-mode 2 \
    --min-seq-id 0.95 \
    --threads 4 \
    --cov-mode 1 \
    -c 0.9 \
    --spaced-kmer-mode 0 \
    --alignment-mode 3 \
    --cluster-reassign 1

mmseqs createtsv genes.db genes.db genes.cluster genes.cluster.tsv
mmseqs version > mmseqs.version
cd ..
```

**Name the clusters and write the catalog** [↗](#name-the-clusters-and-write-the-catalog)

```bash
Rscript solutions/name_clusters.R catalog/genes.cluster.tsv catalog/genes.fna.gz catalog/clusters.tsv
cut -f 2 catalog/clusters.tsv | tail -n +2 > catalog/representatives.txt
seqtk subseq catalog/genes.fna.gz catalog/representatives.txt > catalog/representatives.fna
seqtk subseq catalog/genes.faa.gz catalog/representatives.txt > catalog/representatives.faa
for ext in fna faa; do
  awk 'NR==FNR { if (FNR > 1) name[$2] = $1; next }
       /^>/ { id = substr($1, 2); print ">" name[id]; next } { print }' \
      catalog/clusters.tsv catalog/representatives.$ext > catalog/catalog.$ext
done
head -4 catalog/clusters.tsv
```

:::
