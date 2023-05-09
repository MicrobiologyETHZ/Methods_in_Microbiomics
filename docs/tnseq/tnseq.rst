===================
TN-Seq Protocol
===================
.. _Tim van Opijnen et al.: https://doi.org/10.1038/nmeth.1377
As descibed by `Tim van Opijnen et al.`_, Tn-seq is a method to determine quantitative genetic interactions on a genome-wide scale in microorganisms.

-----------------
Sequence Mapping
-----------------
In the first step of the data processing, we use the reads as they come from the Illumina machine, we identify in the
reads if they carry a Tn5-end; and select only those reads using cutadapt. Then we blast the reads to the transposon
sequence to remove those that are only transposon. The remaining reads are mapped to the combined genome using
bowtie, counted per position and grouped per gene in a perl script.

---------------------
Combined clean Reads
---------------------
In the second step of the data processing, we combine all reads from all samples per library, remove outliers, normalize
and save the files.

-------------
PCA Analysis
-------------
In the third step, we compare all data in a PCA. For this, replace all NaNs by 0, calculate the pca and plot the first
two coefficients.

----------------------
Essentiality Analysis
----------------------
In the third step, we try to define the essential genes at time t=0. We do this by comparing to a randomized virtual
insertion library.

To create the random insertion library, we take the combined genomes and use the actual gene coordinates to outline
starts and stops. We produce a random genome file and choose four times randomly unique values within this range. This
is then summed for the gene positions to create fake random insertion libraries.

To get the t=0 dataset, select the t=0 samples in the normalized files with outliers removed. Replace the NaN by zeros
and estimate whether the mean read number is significantly lower in T0 datasets.

Next, we attribute the group of essentials to their Clusters of Orthologous Genes Category (COG) and compare this
attribution to a random simulation of genes drawn. This we do 10 times, to get an average random COG attribution of
genes from the genome. From there, we calculate the p-value for the difference of the observed to the simulated
COG attributions and plot this as a histogram.

----------------
Paired Analysis
----------------
In the next step, we define the genes with conditional fitness effects at any of the three sampling time points compared
to time t=0. We do this by taking the means across all replicates of the T1-T3 samples, and calculating the ratio to the
mean of the T0 values. We do this for each library individually and for each condition. From this comparison, we also
calculate a p-value and an FDR, from the function mattest and mafdr using BHFDR correction for multiple hypothesis
testing.

As conditional fitness we put a threshold of a ratio change of >2 (‘down’) or <0.5 (‘up’), plus an FDR<0.05, and the
initial ratio to the random insertion library >0.25 (to avoid taking previously allocated essential genes).

Finally, we compare gene lists to find the overlap and unique differences per library and condition. This is used to
plot the Venn diagram.

For defined sets of gene groups we plot the ratio of mean counts at each time point compared to T0 in colour, overlaid
on all ratios; to highlight time trends.

-----------------
Gene Group Plots
-----------------
In this step, we use free text searching and defined gene groups for biological systems such as ‘flagella’ to select the
gene insertion count ratios and plot those. We use both heatmap representation for individual genes, as well as a summed
representation in which all counts for the group is combined and compared.


------------------
COG KEGG Analysis
------------------
In this step of the analysis, we focus analysis on attribution of genes from the selected UP and DOWN gene lists to COG
or KEGG groups, and calculate statistical enrichment or depletion of such group attribution compared to random models.
Finally, we plot the data in a variety of ways.

COG attribution
----------------
We start with the paired comparison data of the ratios T0-replicates/T1-3 replicates and the corresponding p-values and
FDR. We also load the annotation files of the genomes, the manually curated COG list, the gene name list and the gene
number list. For the attribution we use the selected gene lists and the total DOWN and UP files for the two soils.

We now first build a set of random simulations from the length of the selected set. For example, if a DOWN file has 820
items, we simulate random picking for 820 genes from the genome. For each of the randomly picked gene lists, we find
matches to the COG and retrieve the corresponding distribution of COG classes. This is done 10 times. From here we
calculate the means, standard deviation and variation and we compare to the actual COG attributed list in the DOWN file.
Finally, we attribute a confidence interval on the difference to the mean from the random data set by bootstrapping and
a p-value of<0.01, and we calculate a p-value and FDR by mafdr of the random to real list comparison. The results are
saved to a separate Excel for each comparison, to which we also add the individual gene annotations. The whole procedure
and results is also saved in a file.

Polar Plots of COG class distributions
---------------------------------------
All the produced COG class distributions and enrichment factors are now plotted as a polar plot. Essentially, we plot
the mean random relative count of class attribution and the actual relative count, as a log2 transformed bar in a
radial output. We add a star if the FDR value is below 0.05.

KEGG class distributions
-------------------------
Essentially we follow the same idea as for the COG attribution, but now we focus on prediction of metabolic pathways.
Again we make 10 random distributions to compare the observed data to.

Plotting KEGG class enrichment as Volcano plot
-----------------------------------------------
We plot a figure with four panels for each of the exclusive DOWN categories side by side that shows on the x-axis the
log-FC of the observation compared to the mean of the random picked/attributed pathways, and on the y-axis the -log10 of
the calculated p-value for the variation. We also plot a legend with the precise color attribution for the categories.

Plotting metabolic pathways on the KEGG map
--------------------------------------------
Finally, we plot all and exclusive attributed KEGG metabolic pathways on top of the general map in iPATH3, per library
and per soil or liquid category. By gene number or gene name we recover the reaction numbers from the genome-scale model
file, copy that list to the iPATH3 website and save the map as .svg.

