import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

SAMPLES = [
    'aR_LMNB1_25_05_S4',
    'aR_LMNB1_28_06_S10',
    'aR_LMNB1_03_07_S12',
    'nR_LMNB1_21_06_S9',
    'nR_LMNB1_24_05_S1',
    'nR_LMNB1_9_6_S7'
]
SAMPLE_NAMES = [
    'aR25',
    'aR28',
    'aR3',
    'nR21',
    'nR24',
    'nR9'
]

softVSstiff_samples = [
    'Nili_20230219_ExStiff_Ab16048_S13',
    'Nili_20230507_ExStiff_Ab16048_S17',
    'Nili_JAN252023_Exstiff_Ab16048_S19',
    'Nili_20230219_Soft_Ab16048_S11',
    'Nili_20230507_Soft_Ab16048_S15',
    'Nili_JAN252023_Soft_Ab16048_S9',
]
softVSstiff_sample_names = [
    'stiff_S13',
    'stiff_S17',
    'stiff_S19',
    'soft_S11',
    'soft_S15',
    'soft_S9',
]

all_names = [
    'stiff_S13',
    'stiff_S17',
    'stiff_S19',
    'soft_S11',
    'soft_S15',
    'soft_S9',
    'nR21',
    'aR28',
    'aR3',
    'aR25',
    'nR24',
    'nR9'
]



# X:21, Y:22, MT:23

path = 'C:\\Users\\lahat\\Documents\\Uni\\Year4\\Buxboim Lab\\epic2\\2023-24-01\\{}.clip.sort.bed'
chrom_sizes = 'C:\\Users\\lahat\\Documents\\Uni\\Year4\\Buxboim Lab\\chrom.sizes'
anno = 'C:\\Users\\lahat\\Documents\\Uni\\Year4\\Buxboim Lab\\genes\\anno.txt'
window_size = 100000


def create_chr_size_df():
    with open(chrom_sizes, 'r') as f:
        df = pd.read_csv(chrom_sizes, names=['chr', 'size'], delimiter='\t')
    df = df.append(pd.DataFrame({'chr': '0', 'size': 0}, index=[0]))
    df.replace({'X': '21', 'Y': '22', 'MT': '23'}, inplace=True)
    df = df.astype('int64')
    df.set_index('chr', inplace=True)
    for i in range(1, 24):
        df.at[i, 'size'] = df.at[i, 'size'] + df.at[i - 1, 'size']
    # print(df)
    return df


def encode_sample(sample: str, chr_df):
    sample_path = path.format(sample)
    with open(sample_path, 'r') as sample_f:
        sample_df = pd.read_csv(sample_f, names=['chr', 'start', 'end', 'log2fc'], delimiter='\t')
    sample_df.replace({'X': '21', 'Y': '22', 'MT': '23'}, inplace=True)
    encode_start = lambda start, chrom: start + chr_df.at[int(chrom) - 1, 'size']
    encode_end = lambda end, chrom: end + chr_df.at[int(chrom) - 1, 'size']
    sample_df['encode_start'] = pd.Series(map(encode_start, sample_df['start'], sample_df['chr']))
    sample_df['encode_end'] = pd.Series(map(encode_end, sample_df['end'], sample_df['chr']))
    return pd.DataFrame({'start': sample_df['encode_start'], 'end': sample_df['encode_end']}).sort_values(by=['start'])


def sample_lad_windows(sample_encoding, n_windows):
    one_less = window_size - 1
    lad_binary = []
    for i in range(n_windows):
        window_start = i * window_size
        window_end = window_start + one_less
        slice = sample_encoding[(window_start <= sample_encoding['end']) &
                                (sample_encoding['start'] <= window_end)]
        if not slice.empty:
            lad_binary.append(1)
        else:
            lad_binary.append(0)
    return lad_binary


def create_genome_lad_encoding():
    chr_df = create_chr_size_df()
    rn7_genome_size = chr_df.loc[23]['size']
    n_windows = int(rn7_genome_size / window_size)
    encode_df = pd.DataFrame(index=SAMPLES, columns=range(n_windows))
    for sample in SAMPLES:
        sample_encoding = encode_sample(sample, chr_df)
        sample_lad_window_vec = sample_lad_windows(sample_encoding, n_windows)
        encode_df.loc[sample] = sample_lad_window_vec
    encode_df.to_csv('genome_encoding.csv')
    print('genome encoding csv created')


def create_gene_lad_encoding():
    with open(anno) as f:
        all_genes = pd.read_csv(f, sep='\t', usecols=['Gene_name'])
    encode_df = pd.DataFrame(index=SAMPLES, columns=all_genes.values)
    for sample in SAMPLES:
        sample_genes = []
        with open(f'{sample}_genes.txt') as f:
            content = f.read()
            values_list = content.split('\n')
            values_list = [value for value in values_list if value]
            for gene in all_genes.values:
                gene = gene[0]
                if gene in values_list:
                    sample_genes.append(1)
                else:
                    sample_genes.append(0)
        encode_df.loc[sample] = sample_genes
    encode_df.to_csv('gene_encoding.csv')


def create_gene_lad_encoding_soft_stiff():
    with open(anno) as f:
        all_genes = pd.read_csv(f, sep='\t', usecols=['Gene_name'])
    encode_df = pd.DataFrame(index=softVSstiff_samples, columns=all_genes.values)
    pwd = "C:\\Users\\lahat\\Documents\\Uni\\Year4\\Buxboim Lab\\intersects\\softVSstiff-intersect-genes"
    for sample in softVSstiff_samples:
        sample_genes = []
        with open(f'{pwd}\\{sample}_genes.bed_genes.txt') as f:
            content = f.read()
            values_list = content.split('\n')
            values_list = [value for value in values_list if value]
            for gene in all_genes.values:
                gene = gene[0]
                if gene in values_list:
                    sample_genes.append(1)
                else:
                    sample_genes.append(0)
        encode_df.loc[sample] = sample_genes
    encode_df.to_csv('gene_encoding_soft_stiff.csv')


def pcaing(df):
    _pca = PCA()
    xt = _pca.fit_transform(df)
    xt_df = pd.DataFrame(data=xt[:, :2], columns=['pc1', 'pc2'], index=all_names)
    print('Explained variation per principal component: {}'.format(_pca.explained_variance_ratio_))
    plt.figure(figsize=(10, 10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)
    plt.xlabel('PC1', fontsize=20)
    plt.ylabel('PC2', fontsize=20)
    plt.title("PCA of sample LADs genes", fontsize=20)
    colors = ['red', 'indianred', 'tomato', 'green', 'lime', 'yellowgreen', 'orange',
              'blue','darkblue', 'royalblue', 'goldenrod', 'darkorange']
    for sample, color in zip(all_names, colors):
        plt.scatter(xt_df.loc[sample, 'pc1'],
                    xt_df.loc[sample, 'pc2'],
                    c=color,
                    s=300)
    plt.legend(all_names, prop={'size': 15})
    # plt.savefig('genome_pca.eps', format='eps')
    plt.savefig('gene_pca_all_samples.jpg')
    print('gene pca jpg created')


def genome_pca():
    df = pd.read_csv('genome_encoding.csv', index_col=0)
    pcaing(df)


def gene_pca():
    df = pd.read_csv('gene_encoding_soft_stiff.csv', index_col=0)
    pcaing(df)


def gene_pca_soft_stiff_aged_neonate():
    df1 = pd.read_csv('gene_encoding_soft_stiff.csv', index_col=0)
    df2 = pd.read_csv('gene_encoding.csv', index_col=0)
    df3 = pd.concat([df1, df2])
    pcaing(df3)


if __name__ == '__main__':
    # TODO CHANGE THE PCA FIGURE NAME
    # create_genome_lad_encoding()
    # genome_pca()
    # create_gene_lad_encoding()
    # gene_pca()
    # create_gene_lad_encoding_soft_stiff()
    gene_pca_soft_stiff_aged_neonate()