import os
import sys
import pandas as pd
from Bio import SeqIO

import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def GetColorByNormalizedValue(cmap_name, norm_value):
    if norm_value < 0 or norm_value > 1:
        print("ERROR: value " + str(norm_value) + ' does not belong to [0, 1]')
    cmap =  mplt.colormaps[cmap_name] #plt.cm.get_cmap(cmap_name)
    color = cmap(norm_value)
    return mplt.colors.rgb2hex(color[:3])

def ColorByPercentIdentity(pi, config):
    min_pi = config.pi_min
    max_pi = config.pi_max
    fraction = (min(max(pi, min_pi), max_pi) - min_pi) / (max_pi - min_pi)
    if config.reverse:
        fraction = 1 - fraction
    return GetColorByNormalizedValue(config.cmap, fraction)

def ModifyPos(pos, seq_len, strand):
    if strand == '+':
        return pos
    return seq_len - pos + 1

class Config:
    def __init__(self):
        self.pi_min = 85
        self.pi_max = 100
        self.min_align_len = 5000
        self.plot_scale = 1000
        self.cmap = 'Spectral'
        self.reverse = True
        self.upper_triangle = False
        self.linewidth = 1

def main(data_csv, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    config = Config()

    data_df = pd.read_csv(data_csv)
    seq_list = []
    for i in range(len(data_df)):
        seq = [str(r.seq).upper() for r in SeqIO.parse(data_df['Fasta'][i], 'fasta')][0]
        seq_list.append(seq)
    species_names = list(data_df['Label'])
 
    align_dict = dict()
    for i in range(len(species_names)):
        print('aligning ' + species_names[i] + '...')
        #self dot plot
        out_self = os.path.join(output_dir, 'self_' + str(i) + '-' + species_names[i] + '.tsv')
        if not os.path.exists(out_self):
            os.system('lastz ' + data_df['Fasta'][i] + ' ' + data_df['Fasta'][i] + ' --step=20 --notransition --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,id%  --output=' + out_self)
        print('  ' + out_self + ' was computed')
        align_dict[i, i] = out_self
        # pairwise dot plots
        for j in range(i + 1, len(species_names)):
            out_pair = os.path.join(output_dir, 'pair_' + str(i) + '-' + species_names[i] + '_' + str(j) + '-' + species_names[j] + '.tsv')
            if not os.path.exists(out_pair):
                os.system('lastz ' + data_df['Fasta'][i] + ' ' + data_df['Fasta'][j] + ' --step=20 --notransition --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,id%  --output=' + out_pair)
            print('  ' + out_pair + ' was computed')
            align_dict[i, j] = out_pair

    locus_lens = [len(seq) for seq in seq_list]
    for i in range(len(locus_lens)):
        print(data_df['Label'][i], locus_lens[i])

    min_len = min(locus_lens)
    ratios = [round(l / min_len, 2) for l in locus_lens]
    strands = ['+']
    for i in range(1, len(species_names)):
        align_tsv = align_dict[0, i]
        df = pd.read_csv(align_tsv, sep = '\t')
        max_len = max(df['length1'])
        df_max = df.loc[df['length1'] == max_len].reset_index()
        strands.append(df_max['strand2'][0])
    print('Strands', strands)

    df_dict = dict()
    for idx1, idx2 in align_dict:
        df = pd.read_csv(align_dict[idx1, idx2], sep = '\t')
        df = df.loc[(df['length1'] >= config.min_align_len) & (df['length2'] >= config.min_align_len)].reset_index()
        df['id%'] = [float(df['id%'][i][:-1]) for i in range(len(df))]
        df_dict[idx1, idx2] = df

    #### redirecting alignments
    for idx1, idx2 in align_dict:
        len1 = locus_lens[idx1]
        len2 = locus_lens[idx2]
        df = df_dict[idx1, idx2]
        directed_pos_list = []
        for i in range(len(df)):
            pos1 = df['start1'][i], df['end1'][i]
            pos2 = df['start2+'][i], df['end2+'][i]
            if df['strand2'][i] == '-':
                pos2 = df['end2+'][i], df['start2+'][i]
            pos1 = ModifyPos(pos1[0], len1, strands[idx1]) - 1, ModifyPos(pos1[1], len1, strands[idx1]) - 1
            #pos1 = [len1 - pos1[0], len1 - pos1[1]]
            pos2 = ModifyPos(pos2[0], len2, strands[idx2]) - 1, ModifyPos(pos2[1], len2, strands[idx2]) - 1
            directed_pos_list.append([pos1[0], pos1[1], pos2[0], pos2[1]])
        df['start1_dir'] = [p[0] for p in directed_pos_list]
        df['end1_dir'] = [p[1] for p in directed_pos_list]
        df['start2_dir'] = [p[2] for p in directed_pos_list]
        df['end2_dir'] = [p[3] for p in directed_pos_list]

    #### summary alignment stats
    stats_df = {'Label1' : [], 'Label2' : [], 'Idx1' : [], 'Idx2' : [], 'PI' : []}
    for idx1, idx2 in df_dict:
        if idx1 == idx2:
            continue
        df = df_dict[idx1, idx2]
        for i in range(len(df)):
            stats_df['Label1'].append(data_df['Label'][idx1])
            stats_df['Idx1'].append(idx1)
            stats_df['Label2'].append(data_df['Label'][idx2])
            stats_df['Idx2'].append(idx2)
            stats_df['PI'].append(df['id%'][i])
    stats_df = pd.DataFrame(stats_df)
    stats_df.to_csv(os.path.join(output_dir, 'stats_pi.csv'), index = False)

    #### plot setup
    width_ratios = ratios + [0.2]
    if not config.upper_triangle:
        width_ratios = [0.2] + ratios

    fig, axes = plt.subplots(nrows = len(species_names), ncols = len(species_names) + 1, figsize = (20, 20), gridspec_kw={'height_ratios': ratios, 'width_ratios' : width_ratios})

    for i in range(len(data_df)):
        for j in range(len(data_df) + 1):
            plt.xticks([], [])
            plt.yticks([], [])
            axes[i, j].axis("off")

    #### plotting alignments
    for idx1, idx2 in align_dict:
        if config.upper_triangle:
            plt.sca(axes[idx1][idx2])
            axes[idx1][idx2].axis('on')
        else:
            plt.sca(axes[idx2][idx1 + 1])
            axes[idx2][idx1 + 1].axis('on')
        len1 = locus_lens[idx1]
        len2 = locus_lens[idx2]
        df = df_dict[idx1, idx2]
        for i in range(len(df)):
            pos1 = [df['start1_dir'][i], df['end1_dir'][i]]
            pos2 = [df['start2_dir'][i], df['end2_dir'][i]]
            scaled_pos1 = [pos / len1 * config.plot_scale for pos in pos1]
            scaled_pos2 = [pos / len2 * config.plot_scale for pos in pos2]
            pi = df['id%'][i]
            pi_color = ColorByPercentIdentity(pi, config)
            if config.upper_triangle:
                plt.plot([scaled_pos2[0], scaled_pos2[1]], [config.plot_scale - scaled_pos1[0], config.plot_scale - scaled_pos1[1]], color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None')
            else:
                plt.plot([scaled_pos1[0], scaled_pos1[1]], [config.plot_scale - scaled_pos2[0], config.plot_scale - scaled_pos2[1]], color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None')
        plt.xlim(0, config.plot_scale)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])

    gene_col_idx = len(data_df)
    if not config.upper_triangle:
        gene_col_idx = 0
    #### adding V genes
    for idx in range(len(data_df)):
        start_pos = data_df['StartPos'][idx]
        gene_df = pd.read_csv(data_df['GeneTxt'][idx], sep = '\t')
        gene_df['LocusPos'] = [gene_df['Pos'][i] - start_pos for i in range(len(gene_df))]
        locus_len = locus_lens[idx]
        plt.sca(axes[idx][gene_col_idx])
        axes[idx][gene_col_idx].axis('on')
        for i in range(len(gene_df)):
            gene_pos = ModifyPos(gene_df['LocusPos'][i], locus_len, strands[idx])
            scale_pos = config.plot_scale - gene_pos / locus_len * config.plot_scale
            #if config.upper_triangle:
            #    scale_pos = config.plot_scale - scale_pos
            plt.plot([0, 1], [scale_pos, scale_pos], linestyle = '-', marker = 'None', color = 'black')
        plt.xlim(0, 1)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])

    plt.subplots_adjust(hspace=0, wspace = 0)
    plt.savefig(os.path.join(output_dir, '_dotplot.png'), dpi = 300) #, transparent = True)

if __name__ == '__main__':
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    main(input_csv, output_dir)
