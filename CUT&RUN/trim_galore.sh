#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv/trim_galore/%x.%j.out

### Displays job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Job started at `date +"%T %s %d %b %Y"`

### Variables
DIRECTORY=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv
PATH_R1_FASTQ=$1
PATH_R2_FASTQ=$2
NAME=$(basename $1 _L001_R1_001.fastq)
echo $NAME
FULL_PATH="$DIRECTORY/trim_galore/$NAME"

### Modules
module load trim-galore

### Create ouput dir for trimmed fastqs
rm -r $FULL_PATH
mkdir $FULL_PATH

### Adaptor removal with trim galore
trim_galore --illumina --paired --fastqc \
	-o $FULL_PATH \
	$PATH_R1_FASTQ \
	$PATH_R2_FASTQ

### Final time stamp
echo Job finished at `date +"%T %s %d %b %Y"`

