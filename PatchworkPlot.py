import os
import sys
import pandas as pd
from Bio import SeqIO

import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sys.path.append('py')
import visualization_builder as vis
import utils
import data_utils

class Config:
    def __init__(self):
        self.pi_min = 90
        self.pi_max = 95
        self.min_align_len = 10000
        self.plot_scale = 1000
        self.cmap = 'Spectral'
        self.reverse = True
        self.upper_triangle = True
        self.linewidth = 1
        self.transparent = False

def VisualizePlot(plot_utils, ratios, aligned_data, config):
    #### plot setup
    width_ratios = plot_utils.GetWidthRatios(ratios)
    fig, axes = plt.subplots(nrows = plot_utils.NumRows(), ncols = plot_utils.NumColumns(), figsize = (20, 20), gridspec_kw={'height_ratios': ratios, 'width_ratios' : width_ratios})

    #### setting up axes
    for i in range(plot_utils.NumRows()):
        for j in range(plot_utils.NumColumns()):
            plt.xticks([], [])
            plt.yticks([], [])
            axes[i, j].axis("off")

    #### plotting alignments
    for idx1, idx2 in aligned_data.IndexPairIterator():
        plot_utils.SetCurrentAxes(axes, idx1, idx2)
        df = aligned_data.GetAlignmentDF(idx1, idx2)
        len1 = aligned_data.GetLengthByIdx(idx1)
        len2 = aligned_data.GetLengthByIdx(idx2)
        for i in range(len(df)):
            pos1 = [df['start1_dir'][i], df['end1_dir'][i]]
            pos2 = [df['start2_dir'][i], df['end2_dir'][i]]
            scaled_pos1 = [pos / len1 * config.plot_scale for pos in pos1]
            scaled_pos2 = [pos / len2 * config.plot_scale for pos in pos2]
            pi = df['id%'][i]
            pi_color = utils.ColorByPercentIdentity(pi, config)
            x, y = plot_utils.GetLineCoordinates(scaled_pos1[0], scaled_pos1[1], scaled_pos2[0], scaled_pos2[1], config.plot_scale)
            plt.plot(x, y, color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None') 
        plt.xlim(0, config.plot_scale)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])

    #### adding labels
    plot_utils.SetLabels(axes, [aligned_data.GetLabelByIdx(idx) for idx in range(aligned_data.NumSamples())])

    #### adding V genes
    #for idx in range(len(data_df)):
    #    start_pos = data_df['StartPos'][idx]
    #    gene_df = pd.read_csv(data_df['GeneTxt'][idx], sep = '\t')
    #    gene_df['LocusPos'] = [gene_df['Pos'][i] - start_pos for i in range(len(gene_df))]
    #    locus_len = locus_lens[idx]
    #    plt.sca(axes[idx][gene_col_idx])
    #    axes[idx][gene_col_idx].axis('on')
    #    for i in range(len(gene_df)):
    #        gene_pos = utils.ModifyPos(gene_df['LocusPos'][i], locus_len, strands[idx])
    #        scale_pos = config.plot_scale - gene_pos / locus_len * config.plot_scale
    #        plt.plot([0, 1], [scale_pos, scale_pos], linestyle = '-', marker = 'None', color = 'black')
    #    plt.xlim(0, 1)
    #    plt.ylim(0, config.plot_scale)
    #    plt.xticks([], [])
    #    plt.yticks([], [])

    plt.subplots_adjust(hspace = 0, wspace = 0)
    plt.savefig(os.path.join(output_dir, '_dotplot.png'), dpi = 300, transparent = config.transparent)
    plt.savefig(os.path.join(output_dir, '_dotplot.pdf'), dpi = 300)
    plt.clf()

def GetRatios(aligned_data):
    locus_lens = [aligned_data.GetLengthByIdx(i) for i in range(aligned_data.NumSamples())]
    min_len = min(locus_lens)
    ratios = [round(l / min_len, 2) for l in locus_lens]
    return ratios

def main(data_csv, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    config = Config()

    input_data = data_utils.InputData(data_csv)
    pairwise_aligner = data_utils.LastZPairwiseAligner(output_dir)
    aligned_data = data_utils.AlignedData(input_data, pairwise_aligner, output_dir, config)
    aligned_data.ReportSummaryAlignmentStats(os.path.join(output_dir, 'pi_stats.csv')) 

    ratios = GetRatios(aligned_data)
    gene_visualizer = vis.EmptyGeneVisualizer(config, aligned_data.NumSamples())
    plot_visualizer = vis.UpperTriangleUtils(aligned_data.NumSamples(), gene_visualizer)
    VisualizePlot(plot_visualizer, ratios, aligned_data, config)


if __name__ == '__main__':
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    main(input_csv, output_dir)
