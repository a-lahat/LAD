#!/bin/bash
#SBATCH --job-name=picard
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv/picard/%x.%j.out

### Displays job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Job started at `date +"%T %s %d %b %Y"`

### Variables
DIRECTORY=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv
PATH_SORTED_BAM=$1
NAME=$(basename $1 _bowtie2.sorted.bam)
FULL_PATH="$DIRECTORY/picard/$NAME"
OUTPUT_FILE="/${NAME}.modified_bam.bam"
METRICS_FILE="/${NAME}.metrics.txt"

### Create ouput dir for bam after removing duplicates
rm -r $FULL_PATH
mkdir $FULL_PATH

### Modules
module load picard
module load samtools

### Remove duplicate reads with Picard
picard MarkDuplicates \
		I=$PATH_SORTED_BAM \
		O="$FULL_PATH$OUTPUT_FILE" \
		M="$FULL_PATH$METRICS_FILE" \
		REMOVE_DUPLICATES=true

### index bams for epic2
samtools index "$FULL_PATH$OUTPUT_FILE" -o "$FULL_PATH$OUTPUT_FILE.bai"

### Final time stamp
echo Job finished at `date +"%T %s %d %b %Y"`

