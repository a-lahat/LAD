# LAD
## CUT&RUN LAD Profiling Pipeline
For each control and treatment paired-end files of each sample, trim adaptors using **trim-galore**. This is done using a bash script I created that runs the trim-galore module. The script creates an output directory for each control and treatment paired-end files of each sample and writes the output of the trim-galore module to there.

For each control and treatment paired-end files of each sample, map the paired-end reads to the rn7 genome. The genome index is found in the /sci/data/… and the mapping algorithm used is **Bowtie2**. After the mapping, using samtools, the sam file is converted to bam, then the bam is sorted and indexed. 

For each control and treatment mapped, sorted and indexed bam file of each sample, duplicate reads are removed using **picard**. The modified bam is then indexed again using samtools. 

For each sample, find LADs using broad peak-caller **EPIC2**. The EPIC2 tool was installed from github.

For each sample, create bigwig tracks from the bed file output of EPIC2. First I used **bedClip** to clip the peaks by the chrom.sizes file. Then I used bedSort to sort the new bed file. Finally, I used **bedGraphToBigWig** to create the bigwig files. 

To create the cmd for each of these bash scripts I use **runcmds.py** 

After each step in the peak-calling pipeline I check the number of reads in the output to validate that nothing drastic happened.

**Sanity checks:**
I expected the number of reads in the original fastqs and output of trim-galore to have a similar number of reads (checked with wc -l / 4).

I expected the bowtie output to have ~2X as many reads as the trim-galore output because it maps the paired-end reads.

I expected the picard output to have ~80% of the reads as the bowtie output because it removes duplicate reads. (checked with samtools view -c -F 260)

**EPIC2 parameters** - choosing chromosomes (canonical?), optimizing parameters - bin width, gap size

**LAD PCA:**
LAD genes PCA - use bedtools intersect to intersect epic2 output bed with rn7 annotation file to find which genes fall in LADs. This created a binary matrix of columns=genes, rows=samples, that describes which genes appear in LADs in which samples. PCA on this matrix shows us samples that are close together. Neonate samples are close together, but aged samples are far from neonates and from each other. This can be explained by - each cell goes through different things in life so it has different LADs. The question is how much is common between the aged sample LADs.

LAD genome PCA - same as above but divide genome to bins of size X and for each bin mark 1 if this bin has a part of a LAD in it, 0 if not. 

**Defining new LADs** - 10% overlap between LADs is considered the same LAD, so rewrite the bed files of epic2 output accordingly.

Same steps for soft VS stiff samples.

**CUT&RUN data gives us** - list of LADs per sample, list of genes in LADs per sample, percentage of genome covered by LADs per sample, PCA of sample.

## Single-cell RNA Sequencing Analysis Pipeline
The scRNA-seq data is a matrix where the columns are cells, the rows are genes, and each cell in the table is the amount of expression of that gene in that cell.

Using the annotation file for rn7, I switch the gene stable id with the gene name if it’s known. Then I remove all the genes that are MT- genes.

For the scRNA-seq analysis I used Scanpy, a library for single cell analysis in python based a lot on Seurat. 

Preprocessing - keep genes that are expressed in at least 3 cells, keep cells that express at least 200 genes.

Normalization - normalizes total counts per cell to 10,000, logarithmizes the data matrix by computing X = ln(X + 1), regresses out (mostly) unwanted sources of variation using simple linear regression, and scales data to unit variance and zero mean and truncate values to 10.

Pca - Uses the implementation of scikit-learn

Umap - First, compute a neighborhood graph of observations, then embed the neighborhood graph using UMAP and cluster cells using the Leiden algorithm.

UMAP can be colored by - sample types (aOPC, nOPC, soft, stiff), sample labels (each sample by itself), leiden → cell types by gene markers, specific gene expression intensity (black&white graphs) to show different things. 

The next step was to use the final matrix from the sc analysis to find the mean expression of each gene in each defined cluster. Then, using a list of rn7 neuronal gene markers, I created a matrix of columns=clusters, rows=gene markers. For each cluster I identified which cell type it is by the gene marker(s) that were most expressed.

Use the scanpy rank_gene_groups method to get the log2fc and pvalue for each gene between the groups nOPC VS aOPC.

Optimizing the **intensity threshold** by creating a histogram of the intensities and choosing the second peak in the bimodule graph. (sum_hist.png)

Optimize the **thresholds** for **pvalue** and **log2foldchange** by iterating through multiple thresholds and creating a scatter plot of the ratio of aged-exclusive LAD genes that are down-regulated VS  the ratio of neonate-exclusive LAD genes that are up-regulated (which_log2.png). Choosing the threshold with the highest ratios.

Define up/down-regulated genes as genes that pass **all** these thresholds.

Find remyelination related genes in the up/down-regulated gene groups.

See how many of the genes in the **aged**-exclusive LAD group are in the **down**-regulated group, how many of the genes in the **neonate**-exclusive LAD group are in the **up**-regulated group.

Tweak thresholds to make ratios significant.

**SCRNA-SEQ data gives us** - up/down-regulated genes, gene expression differences in cell types/sample types.
