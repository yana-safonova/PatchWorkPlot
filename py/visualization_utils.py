import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colorbar import ColorbarBase

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

    def VisualizePairwiseGenes(self, axes, idx1, idx2):
        axes[1, 0].axis("off")
        axes[0, 1].axis("off")

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
            plt.xlim(0, 1)
            plt.ylim(0, self.config.plot_scale)
            plt.xticks([], [])
            plt.yticks([], [])

    def VisualizePairwiseGenes(self, axes, idx1, idx2):
        plt.sca(axes[1, 0])
        gene_df2 = self.aligned_data.GetGeneTableByIdx(idx2)
        locus_len2 = self.aligned_data.GetLengthByIdx(idx2)
        strand2 = self.aligned_data.GetStrandByIdx(idx2)
        for i in range(len(gene_df2)):
            pos_list = utils.ModifyPos(gene_df2['Start'][i], locus_len2, strand2), utils.ModifyPos(gene_df2['End'][i], locus_len2, strand2)
            scaled_pos_list = [gene_pos / locus_len2 * self.config.plot_scale for gene_pos in pos_list]
            width = max(2, abs(scaled_pos_list[0] - scaled_pos_list[1]))
            rect = patches.Rectangle((min(scaled_pos_list), 0), width, 1, linewidth=0, edgecolor='r', facecolor=gene_df2['Color'][i])
            plt.gca().add_patch(rect)
        plt.xlim(0, self.config.plot_scale)
        plt.ylim(0, 1)
        plt.xticks([], [])
        plt.yticks([], [])
        ####
        if idx1 == idx2:
            axes[0, 1].axis("off")
            return
        plt.sca(axes[0, 1])
        gene_df1 = self.aligned_data.GetGeneTableByIdx(idx1)
        locus_len1 = self.aligned_data.GetLengthByIdx(idx1)
        strand1 = self.aligned_data.GetStrandByIdx(idx1)
        for i in range(len(gene_df1)):
            pos_list = utils.ModifyPos(gene_df1['Start'][i], locus_len1, strand1), utils.ModifyPos(gene_df1['End'][i], locus_len1, strand1)
            scaled_pos_list = [(self.config.plot_scale - gene_pos / locus_len1 * self.config.plot_scale) for gene_pos in pos_list]
            width = max(2, abs(scaled_pos_list[0] - scaled_pos_list[1]))
            rect = patches.Rectangle((0, min(scaled_pos_list)), 1, width, linewidth=0, edgecolor='r', facecolor=gene_df1['Color'][i])
            plt.gca().add_patch(rect)
        plt.xlim(0, 1)
        plt.ylim(0, self.config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])


class ColorUtils:
    def __init__(self, config):
        self.config = config

    def GetColor(self, pi):
        if self.config.color != '':
            return self.config.color
        return utils.ColorByPercentIdentity(self.config.cmap, pi, self.config.pi_min, self.config.pi_max, self.config.cmap_reverse)

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

    def VisualizePairwiseGenes(self, axes, idx1, idx2):
        self.gene_vis_utils.VisualizePairwiseGenes(axes, idx1, idx2)

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

    def VisualizePairwiseGenes(self, axes, idx1, idx2):
        self.gene_vis_utils.VisualizePairwiseGenes(axes, idx1, idx2)


def GetRatios(aligned_data):
    locus_lens = [aligned_data.GetLengthByIdx(i) for i in range(aligned_data.NumSamples())]
    min_len = min(locus_lens)
    ratios = [round(l / min_len, 2) for l in locus_lens]
    return ratios

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
    color_utils = ColorUtils(config)
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
            pi_color = color_utils.GetColor(pi)
            x, y = plot_utils.GetLineCoordinates(scaled_pos1[0], scaled_pos1[1], scaled_pos2[0], scaled_pos2[1], config.plot_scale)
            plt.plot(x, y, color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None')
            #### add breakpoints
            if config.show_breakpoints:
                if min(abs(pos2[1] - pos2[0]), abs(pos1[1] - pos1[0])) <= config.bp_min_len:
                    continue
                plt.plot([0, config.plot_scale], [y[0], y[0]], color = config.bp_color, linestyle = '-', marker = 'None', linewidth = config.bp_linewidth)
                plt.plot([0, config.plot_scale], [y[1], y[1]], color = config.bp_color, linestyle = '-', marker = 'None', linewidth = config.bp_linewidth)
                plt.plot([x[0], x[0]], [0, config.plot_scale], color = config.bp_color, linestyle = '-', marker = 'None', linewidth = config.bp_linewidth)
                plt.plot([x[1], x[1]], [0, config.plot_scale], color = config.bp_color, linestyle = '-', marker = 'None', linewidth = config.bp_linewidth)
        plt.xlim(0, config.plot_scale)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])

    #### add legend
    if not config.hide_legend:
        if config.upper_triangle:
            cax = fig.add_axes([0.64, 0.09, 0.258, 0.007])
        else:
            cax = fig.add_axes([0.125, 0.09, 0.258, 0.007])

        if config.cmap_reverse:
            cmap = plt.get_cmap(config.cmap).reversed()
        else:
            cmap = plt.get_cmap(config.cmap)

        ColorbarBase(
            cax,
            cmap=cmap,
            norm=plt.Normalize(config.pi_min, config.pi_max),
            orientation='horizontal',
            ticks=[config.pi_min, config.pi_max],
            label='% identity'
        )

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
    plt.savefig(os.path.join(config.output_dir, 'patchworkplot.png'), dpi = 300, bbox_inches='tight',transparent = config.transparent)
    plt.savefig(os.path.join(config.output_dir, 'patchworkplot.pdf'), dpi = 300, bbox_inches='tight',)
    plt.clf()

def GetFigureSizes(len1, len2):
    max_len = 6
    min_len = 6 * min(len1, len2) / max(len1, len2)
    if len1 == max(len1, len2):
        return min_len, max_len
    return max_len, min_len

def PlotPairwiseAlignments(plot_utils, aligned_data, config):
    color_utils = ColorUtils(config)
    ratio = 10
    for idx1 in range(aligned_data.NumSamples()):
        for idx2 in range(idx1, aligned_data.NumSamples()):
            df = aligned_data.GetAlignmentDF(idx1, idx2)
            locus_len1 = aligned_data.GetLengthByIdx(idx1)
            locus_len2 = aligned_data.GetLengthByIdx(idx2)
            s1, s2 = GetFigureSizes(locus_len1, locus_len2)
            ratio2 = s1 * (ratio + 1) / s2
            fig, axes = plt.subplots(2, 2, figsize = (s1, s2), gridspec_kw={'height_ratios': [ratio, 1], 'width_ratios' : [(ratio2 + 1), 1]})
            plt.sca(axes[0, 0])
            for i in range(len(df)):
                pos1 = [df['start1_dir'][i], df['end1_dir'][i]]
                pos2 = [df['start2_dir'][i], df['end2_dir'][i]]
                scaled_pos1 = [pos / locus_len1 * config.plot_scale for pos in pos1]
                scaled_pos2 = [pos / locus_len2 * config.plot_scale for pos in pos2]
                pi = df['id%'][i]
                pi_color = color_utils.GetColor(pi)
                x, y = plot_utils.GetLineCoordinates(scaled_pos1[0], scaled_pos1[1], scaled_pos2[0], scaled_pos2[1], config.plot_scale)
                plt.plot(x, y, color = pi_color, linewidth=config.linewidth, linestyle = '-', marker = 'None')
            plt.xlim(0, config.plot_scale)
            plt.ylim(0, config.plot_scale)
            plt.xticks([], [])
            plt.yticks([], [])
            plot_utils.VisualizePairwiseGenes(axes, idx1, idx2)
            axes[1, 1].axis("off")
            plt.subplots_adjust(hspace = 0, wspace = 0)
            plt.savefig(os.path.join(config.pairwise_plot_dir, str(idx1) + '-' + aligned_data.GetSampleNameByIdx(idx1) + '_' + str(idx2) + '-' + aligned_data.GetSampleNameByIdx(idx2) + '.png'), dpi = 300)
            plt.clf()
            plt.close()
