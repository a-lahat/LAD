#!/bin/bash
#SBATCH --job-name=bowtie2_paired
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv/bowtie2/%x.%j.out

### Displays job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Job started at `date +"%T %s %d %b %Y"`

### Variables
DIRECTORY=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv
MRATBN7_GENOME_INDEX=/sci/data/reference_data/Rattus_norvegicus/Ensembl/mRatBN7.2/Sequence/Bowtie2Index/genome
PATH_R1_FASTQ=$1
PATH_R2_FASTQ=$2
NAME=$(basename $1 _L001_R1_001_val_1.fq)
FULL_PATH="$DIRECTORY/bowtie2/$NAME"
SAM="/${NAME}_bowtie2.sam"
STDERR="/${NAME}_bowtie2.stderr"
BAM="/${NAME}_bowtie2.bam"
SORTED_BAM="/${NAME}_bowtie2.sorted.bam"

### Create ouput dir for mapped files
rm -r $FULL_PATH
mkdir $FULL_PATH

### Map reads with bowtie2
bowtie2 --local --very-sensitive-local --phred33 --no-mixed --no-discordant --no-unal --threads 16 \
	-x $MRATBN7_GENOME_INDEX \
	-I 10 -X 700 \
	-1 $PATH_R1_FASTQ \
	-2 $PATH_R2_FASTQ \
	-S "$FULL_PATH$SAM" \
	2>"$FULL_PATH$STDERR"

### Convert SAM to BAM
samtools view -Sb \
	"$FULL_PATH$SAM" \
	-o "$FULL_PATH$BAM"

### Sort BAM
samtools sort \
       "$FULL_PATH$BAM"	\
	-o "$FULL_PATH$SORTED_BAM"

### Indexes BAM
samtools index -b \
	"$FULL_PATH$SORTED_BAM"

### Final time stamp
echo Job finished at `date +"%T %s %d %b %Y"`

