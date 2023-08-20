#!/bin/bash
#SBATCH --job-name=log2bw
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv/log2bw/%x.%j.out

### Displays job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Job started at `date +"%T %s %d %b %Y"`

### Variables
DIRECTORY=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv
PATH_BED=$1
NAME=$(basename $1 .epic2.bed)
FULL_PATH="$DIRECTORY/log2bw/$NAME"
CHRSIZE="$DIRECTORY/chrom.sizes"

### Create ouput dir for bam after removing duplicates
rm -r $FULL_PATH
mkdir $FULL_PATH

### Modules
module load bedtools

cut -f 1,2,3,10 $PATH_BED >> "$FULL_PATH/${NAME}.log2.bed"

### clip bed
"$DIRECTORY/bedClip" "$FULL_PATH/${NAME}.log2.bed" "${CHRSIZE}" "$FULL_PATH/${NAME}.clip.bed"

### sort clip
"$DIRECTORY/bedSort" "$FULL_PATH/${NAME}.clip.bed" "$FULL_PATH/${NAME}.clip.sort.bed"

### bed to bigwig
"$DIRECTORY/bedGraphToBigWig" "$FULL_PATH/${NAME}.clip.sort.bed" "${CHRSIZE}" "$FULL_PATH/${NAME}.bw"
