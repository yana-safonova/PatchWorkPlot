import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import utils

class EmptyGeneVisualizer:
    def __init__(self, config, aligned_data):
        self.config = config

    def ShowGenesFlag(self):
        return False

    def GetColumnWidth(self):
        return []

    def NumGeneColumns(self):
        return 0

    def VisualizeGenes(self, axes, gene_column_idx):
        return 

class SimpleGeneVisualizer:
    def __init__(self, config, aligned_data):
        self.config = config
        self.aligned_data = aligned_data

    def ShowGenesFlag(self):
        return True

    def GetColumnWidth(self):
        return [0.2]

    def NumGeneColumns(self):
        return 1

    def VisualizeGenes(self, axes, gene_col_idx):
        for idx in range(self.aligned_data.NumSamples()):
            plt.sca(axes[idx][gene_col_idx])
            axes[idx][gene_col_idx].axis('on')
            gene_df = self.aligned_data.GetGeneTableByIdx(idx)
            locus_len = self.aligned_data.GetLengthByIdx(idx)
            strand = self.aligned_data.GetStrandByIdx(idx)
            for i in range(len(gene_df)):
                pos_list = utils.ModifyPos(gene_df['Start'][i], locus_len, strand), utils.ModifyPos(gene_df['End'][i], locus_len, strand)
                scaled_pos_list = [(self.config.plot_scale - gene_pos / locus_len * self.config.plot_scale) for gene_pos in pos_list]
                width = max(2, abs(scaled_pos_list[0] - scaled_pos_list[1]))
                rect = patches.Rectangle((0, min(scaled_pos_list)), 1, width, linewidth=0, edgecolor='r', facecolor=gene_df['Color'][i])
                plt.gca().add_patch(rect)
                #plt.plot([0, 1], [scaled_pos_list[0], scaled_pos_list[1]], linestyle = '-', marker = 'None', color = gene_df['Color'][i])
            plt.xlim(0, 1)
            plt.ylim(0, self.config.plot_scale)
            plt.xticks([], [])
            plt.yticks([], [])            


class UpperTriangleUtils:
    def __init__(self, aligned_data, gene_vis_utils):
        self.aligned_data = aligned_data
        self.num_samples = aligned_data.NumSamples()
        self.gene_vis_utils = gene_vis_utils

    def SetCurrentAxes(self, axes, idx1, idx2):
        plt.sca(axes[idx1][idx2])
        axes[idx1][idx2].axis('on')

    def GetGeneColumnIndex(self):
        return self.num_samples        

    def GetWidthRatios(self, ratios):
        return ratios + self.gene_vis_utils.GetColumnWidth()

    def NumColumns(self):
        return self.num_samples + self.gene_vis_utils.NumGeneColumns()

    def NumRows(self):
        return self.num_samples

    def GetLineCoordinates(self, x1, x2, y1, y2, scale):
        return [y1, y2], [scale - x1, scale - x2]

    def SetLabels(self, axes, sample_labels):
        for idx, label in enumerate(sample_labels):
            plt.sca(axes[0, idx])
            plt.title(label)

    def VisualizeGenes(self, axes):
        gene_col_idx = self.GetGeneColumnIndex()
        self.gene_vis_utils.VisualizeGenes(axes, gene_col_idx)

class LowerTriangleUtils:
    def __init__(self, aligned_data, gene_vis_utils):
        self.aligned_data = aligned_data
        self.num_samples = aligned_data.NumSamples()
        self.gene_vis_utils = gene_vis_utils
        self.col_shift = int(self.gene_vis_utils.ShowGenesFlag())

    def SetCurrentAxes(self, axes, idx1, idx2):
        plt.sca(axes[idx2][idx1 + self.col_shift])
        axes[idx2][idx1 + self.col_shift].axis('on')

    def GetGeneColumnIndex(self):
        return 0

    def GetWidthRatios(self, ratios):
        return self.gene_vis_utils.GetColumnWidth() + ratios

    def NumColumns(self):
        return self.num_samples + self.gene_vis_utils.NumGeneColumns()

    def NumRows(self):
        return self.num_samples

    def GetLineCoordinates(self, x1, x2, y1, y2, scale):
        return [x1, x2], [scale - y1, scale - y2] 

    def SetLabels(self, axes, sample_labels):
        for idx, label in enumerate(sample_labels):
            plt.sca(axes[self.num_samples - 1, self.col_shift + idx])
            plt.xlabel(label)

    def VisualizeGenes(self, axes):
        gene_col_idx = self.GetGeneColumnIndex()
        self.gene_vis_utils.VisualizeGenes(axes, gene_col_idx)


def GetRatios(aligned_data):
    locus_lens = [aligned_data.GetLengthByIdx(i) for i in range(aligned_data.NumSamples())]
    min_len = min(locus_lens)
    ratios = [round(l / min_len, 2) for l in locus_lens]
    return ratios

def PlotSelfAlignments(plot_utils, aligned_data, config):
    for idx in range(aligned_data.NumSamples()):
        df = aligned_data.GetAlignmentDF(idx, idx)
        locus_len = aligned_data.GetLengthByIdx(idx)
        plt.figure(figsize = (6, 6))
        for i in range(len(df)):
            pos1 = [df['start1_dir'][i], df['end1_dir'][i]]
            pos2 = [df['start2_dir'][i], df['end2_dir'][i]]
            scaled_pos1 = [pos / locus_len * config.plot_scale for pos in pos1]
            scaled_pos2 = [pos / locus_len * config.plot_scale for pos in pos2]
            pi = df['id%'][i]
            pi_color = 'black' #utils.ColorByPercentIdentity(config.cmap, pi, config.pi_min, config.pi_max, config.cmap_reverse)
            x, y = plot_utils.GetLineCoordinates(scaled_pos1[0], scaled_pos1[1], scaled_pos2[0], scaled_pos2[1], config.plot_scale)
            plt.plot(x, y, color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None')
        plt.xlim(0, config.plot_scale)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.savefig(os.path.join(config.output_dir, 'selfdotplot_' + aligned_data.GetSampleNameByIdx(idx) + '.png'), dpi = 300)
        plt.clf()

def VisualizePlot(plot_utils, aligned_data, config):
    #### get ratios
    ratios = GetRatios(aligned_data)
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
            pi_color = utils.ColorByPercentIdentity(config.cmap, pi, config.pi_min, config.pi_max, config.cmap_reverse)
            x, y = plot_utils.GetLineCoordinates(scaled_pos1[0], scaled_pos1[1], scaled_pos2[0], scaled_pos2[1], config.plot_scale)
            plt.plot(x, y, color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None')
        plt.xlim(0, config.plot_scale)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])

    #### adding labels and gene positions
    labels = []
    for idx in range(aligned_data.NumSamples()):
        label = ''
        if not pd.isnull(aligned_data.GetLabelByIdx(idx)):
            label = aligned_data.GetLabelByIdx(idx)
        labels.append(label)
    plot_utils.SetLabels(axes, labels)
    plot_utils.VisualizeGenes(axes)

    #### output plot as .PNG, .PDF
    plt.subplots_adjust(hspace = 0, wspace = 0)
    plt.savefig(os.path.join(config.output_dir, 'patchworkplot.png'), dpi = 300, transparent = config.transparent)
    plt.savefig(os.path.join(config.output_dir, 'patchworkplot.pdf'), dpi = 300)
    plt.clf()
