#!/bin/bash
#SBATCH --job-name=bamtobigwig
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv/bamtobigwig/%x.%j.out

### Displays job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Job started at `date +"%T %s %d %b %Y"`

### Variables
DIRECTORY=/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv
PATH_BAM=$1
NAME=$(basename $1 .modified_bam.bam)
FULL_PATH="$DIRECTORY/bamtobigwig/$NAME"
OUTPUT_FILE="/${NAME}.bedGraph"

### Create ouput dir for bam after removing duplicates
rm -r $FULL_PATH
mkdir $FULL_PATH

### Modules
module load bedtools

### bam to bedGraph
bedtools genomecov -pc -ibam $PATH_BAM -bg >> "$FULL_PATH$OUTPUT_FILE"

### sort bedGraph
sort -k1,1 -k2,2n -k3,3n -s "$FULL_PATH$OUTPUT_FILE" > "$FULL_PATH$OUTPUT_FILE.tmp"

bedtools merge -i "$FULL_PATH$OUTPUT_FILE.tmp" -c 4 -d 0 -o max > "$FULL_PATH$OUTPUT_FILE.out"

sort -k1,1 -k2,2n -k3,3n -s "$FULL_PATH$OUTPUT_FILE.out" > "$FULL_PATH$OUTPUT_FILE.sorted"

### bedGraph to bigwig
"$DIRECTORY/bedGraphToBigWig" "$FULL_PATH$OUTPUT_FILE.sorted" /sci/data/reference_data/Rattus_norvegicus/Ensembl/mRatBN7.2/Sequence/WholeGenomeFasta/sizes.txt "$FULL_PATH/${NAME}.bw"

### Final time stamp
echo Job finished at `date +"%T %s %d %b %Y"`

