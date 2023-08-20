GENOME = '/sci/data/reference_data/Rattus_norvegicus/Ensembl/mRatBN7.2/Sequence/Bowtie2Index/genome'
GENOME_FA = '/sci/data/reference_data/Rattus_norvegicus/Ensembl/mRatBN7.2/Sequence/WholeGenomeFasta/genome.fa'

CHR_SIZE = '/sci/data/reference_data/Rattus_norvegicus/Ensembl/mRatBN7.2/Sequence/WholeGenomeFasta/sizes.txt'
SOME_CHR_SIZES = '/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/ExperimentCut_Run/20220727_NGSData/Nili-358699458/Nili_Patrick_Haneen-589917330/aviv/sizes.txt'

LAB_PATH = '/cs/labs/buxi/Nili.Avidan/DamID/DamID_Experiments/' \
           'ExperimentCut_Run/20220727_NGSData/Nili-358699458/' \
           'Nili_Patrick_Haneen-589917330/aviv'

FASTQS = 'fastqs'
TRIMGALORE = 'trim_galore'
BOWTIE2 = 'bowtie2'
PICARD = 'picard'

aR_samples = [
    ['aR_IgG_25_05_S5',
     'aR_LMNB1_25_05_S4'],
    ['aR_IgG_28_06_S11',
     'aR_LMNB1_28_06_S10'],
    ['aR_C42D8_03_07_S13',
     'aR_LMNB1_03_07_S12']
]

nR_samples = [
    ['nR_IgG_12_07_S15',
     'nR_LMNB1_12_07_S14'],
    ['nR_IgG_21_06_S6',
     'nR_LMNB1_21_06_S9'],
    ['nR_IgG_24_05_S3',
     'nR_LMNB1_24_05_S1'],
    ['nR_IgG_9_6_S8',
     'nR_LMNB1_9_6_S7']
]

softVSstiff_samples = [
    '20230219_ExStiff_Ab16048_S13',
    '20230219_ExStiff_IgG_S14',
    '20230507_ExStiff_Ab16048_S17',
    '20230507_ExStiff_IgG_S18',
    'JAN252023_Exstiff_Ab16048_S19',
    'JAN252023_Exstiff_IgG_S20',
    '20230219_Soft_Ab16048_S11',
    '20230219_Soft_IgG_S12',
    '20230507_Soft_Ab16048_S15',
    '20230507_Soft_IgG_S16',
    'JAN252023_Soft_Ab16048_S9',
    'JAN252023_Soft_IgG_S10'
]

softVSstiff_samples_paired = [
    ['20230219_ExStiff_Ab16048_S13',
     '20230219_ExStiff_IgG_S14', ],
    ['20230507_ExStiff_Ab16048_S17',
     '20230507_ExStiff_IgG_S18', ],
    ['JAN252023_Exstiff_Ab16048_S19',
     'JAN252023_Exstiff_IgG_S20', ],
    ['20230219_Soft_Ab16048_S11',
     '20230219_Soft_IgG_S12', ],
    ['20230507_Soft_Ab16048_S15',
     '20230507_Soft_IgG_S16', ],
    ['JAN252023_Soft_Ab16048_S9',
     'JAN252023_Soft_IgG_S10'],
]

bam_ext = '.modified_bam.bam'


# 'nR_IgG_12_07_S15/nR_IgG_12_07_S15.modified_bam.bam'


def run_trim_galore():
    with open('trim_galore_commands.txt', 'w') as f:
        for sample in softVSstiff_samples:
            R1 = f'{LAB_PATH}/fastqs/Nili_{sample}_L001_R1_001.fastq'
            R2 = f'{LAB_PATH}/fastqs/Nili_{sample}_L001_R2_001.fastq'
            f.write(f'sbatch scripts/trim_galore.sh {R1} {R2}\n')


def run_bowtie2():
    with open('bowtie2_commands.txt', 'w') as f:
        # for samples in [nR_samples, aR_samples]:
        #     for sample in samples:
        for sample in softVSstiff_samples_paired:
            # control_name = sample[0]
            # treatment_name = sample[1]
            control_name = f'Nili_{sample[1]}'
            treatment_name = f'Nili_{sample[0]}'
            control_R1 = f'{LAB_PATH}/{TRIMGALORE}/{control_name}/{control_name}_L001_R1_001_val_1.fq'
            control_R2 = f'{LAB_PATH}/{TRIMGALORE}/{control_name}/{control_name}_L001_R2_001_val_2.fq'
            treatment_R1 = f'{LAB_PATH}/{TRIMGALORE}/{treatment_name}/{treatment_name}_L001_R1_001_val_1.fq'
            treatment_R2 = f'{LAB_PATH}/{TRIMGALORE}/{treatment_name}/{treatment_name}_L001_R2_001_val_2.fq'
            f.write(f'sbatch scripts/bowtie2_paired.sh {control_R1} {control_R2}\n')
            f.write(f'sbatch scripts/bowtie2_paired.sh {treatment_R1} {treatment_R2}\n')


def check_bowtie2_reads():
    with open('bowtie2_reads.txt', 'w') as f:
        # for samples in [nR_samples, aR_samples]:
        #     for sample in samples:
        for sample in softVSstiff_samples_paired:
            # control_name = sample[0]
            # treatment_name = sample[1]
            control_name = f'Nili_{sample[1]}'
            treatment_name = f'Nili_{sample[0]}'
            control = f'{LAB_PATH}/{BOWTIE2}/{control_name}/{control_name}_bowtie2.sorted.bam'
            treatment = f'{LAB_PATH}/{BOWTIE2}/{treatment_name}/{treatment_name}_bowtie2.sorted.bam'
            f.write(f'samtools view -c -F 260 {control}\n')
            f.write(f'samtools view -c -F 260 {treatment}\n')


def run_picard():
    with open('picard_commands.txt', 'w') as f:
        # for samples in [nR_samples, aR_samples]:
        #     for sample in samples:
        for sample in softVSstiff_samples_paired:
            # control_name = sample[0]
            # treatment_name = sample[1]
            control_name = f'Nili_{sample[1]}'
            treatment_name = f'Nili_{sample[0]}'
            control = f'{LAB_PATH}/{BOWTIE2}/{control_name}/{control_name}_bowtie2.sorted.bam'
            treatment = f'{LAB_PATH}/{BOWTIE2}/{treatment_name}/{treatment_name}_bowtie2.sorted.bam'
            f.write(f'sbatch scripts/picard.sh {control}\n')
            f.write(f'sbatch scripts/picard.sh {treatment}\n')


def check_picard_reads():
    with open('picard_reads.txt', 'w') as f:
        # for samples in [nR_samples, aR_samples]:
        #     for sample in samples:
        for sample in softVSstiff_samples_paired:
            # control_name = sample[0]
            # treatment_name = sample[1]
            control_name = f'Nili_{sample[1]}'
            treatment_name = f'Nili_{sample[0]}'
            control = f'{LAB_PATH}/{PICARD}/{control_name}/{control_name}.modified_bam.bam'
            treatment = f'{LAB_PATH}/{PICARD}/{treatment_name}/{treatment_name}.modified_bam.bam'
            f.write(f'samtools view -c -F 260 {control}\n')
            f.write(f'samtools view -c -F 260 {treatment}\n')


def run_bamtobigwig():
    with open('bamtobigwig_commands.txt', 'w') as f:
        for samples in [nR_samples, aR_samples]:
            for sample in samples:
                control_name = sample[0]
                treatment_name = sample[1]
                control = f'{LAB_PATH}/{PICARD}/{control_name}/{control_name}.modified_bam.bam'
                treatment = f'{LAB_PATH}/{PICARD}/{treatment_name}/{treatment_name}.modified_bam.bam'
                f.write(f'sbatch scripts/bamtobigwig.sh {control}\n')
                f.write(f'sbatch scripts/bamtobigwig.sh {treatment}\n')


def run_epic2bw():
    with open('epic2bw_commands.txt', 'w') as f:
        for samples in [nR_samples, aR_samples]:
            for sample in samples:
                control_name = sample[0]
                treatment_name = sample[1]
                control = f'{LAB_PATH}/{PICARD}/{control_name}/{control_name}{bam_ext}'
                treatment = f'{LAB_PATH}/{PICARD}/{treatment_name}/{treatment_name}{bam_ext}'
                epic_cmd_params = f'--treatment {treatment} ' \
                                  f'--control {control} ' \
                                  f'--chromsizes {SOME_CHR_SIZES} ' \
                                  f'--egf 0.74 ' \
                                  f'--bin-size 8000 ' \
                                  f'--gaps-allowed 6 ' \
                                  f'--guess-bampe ' \
                                  f'--log2fc-bigwig log2fc_{treatment_name}.bw'
                print(f'Control: {control_name}')
                print(f'Treatment: {treatment_name}')
                f.write(f'epic2-bw {epic_cmd_params}\n')


def run_epic2():
    with open('epic2_commands.txt', 'w') as f:
        # for samples in [nR_samples, aR_samples]:
        #     for sample in samples:
        for sample in softVSstiff_samples_paired:
            # control_name = sample[0]
            # treatment_name = sample[1]
            control_name = f'Nili_{sample[1]}'
            treatment_name = f'Nili_{sample[0]}'
            control = f'{LAB_PATH}/{PICARD}/{control_name}/{control_name}{bam_ext}'
            treatment = f'{LAB_PATH}/{PICARD}/{treatment_name}/{treatment_name}{bam_ext}'
            epic_cmd_params = f'--treatment {treatment} ' \
                              f'--control {control} ' \
                              f'--chromsizes {SOME_CHR_SIZES} ' \
                              f'-egf 0.992 ' \
                              f'--bin-size 8000 ' \
                              f'--gaps-allowed 6 ' \
                              f'--guess-bampe ' \
                              f'-o {treatment_name}.epic2.bed'
            print(f'Control: {control_name}')
            print(f'Treatment: {treatment_name}')
            f.write(f'epic2 {epic_cmd_params}\n')
            # print(f'cmd: epic2 {epic_cmd_params}')
            # runcmd.run(['epic2', epic_cmd_params])


def run_log2bw():
    with open('log2bw_commands.txt', 'w') as f:
        for sample in softVSstiff_samples_paired:
            treatment_name = f'Nili_{sample[0]}'
            treatment = f'{LAB_PATH}/epic2/{treatment_name}.epic2.bed'
            f.write(f'sbatch scripts/log2bedtobw.sh {treatment}\n')


if __name__ == '__main__':
    # check_bowtie2_reads()
    # check_picard_reads()
    # run_epic2()
    # run_trim_galore()
    # run_bowtie2()
    # run_picard()
    run_log2bw()
