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
