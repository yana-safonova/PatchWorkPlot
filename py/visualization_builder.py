import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import utils

class EmptyGeneVisualizer:
    def __init__(self, config, num_samples):
        self.config = config
        self.num_samples = num_samples

    def ShowGenesFlag(self):
        return False

    def GetColumnWidth(self):
        return []

    def NumGeneColumns(self):
        return 0

    def VisualizeGenes(self, gene_column_idx):
        return 

class RegularGeneVisualizer:
    def __init__(self, config, num_samples):
        self.config = config
        self.num_samples = num_samples

    def ShowGenesFlag(self):
        return True

    def GetColumnWidth(self):
        return [0.2]

    def NumGeneColumns(self):
        return 1

    def VisualizeGenes(self, data_df, gene_column_idx):
        for idx in range(self.num_samples):
            start_pos = data_df['StartPos'][idx]
            gene_df = pd.read_csv(data_df['GeneTxt'][idx], sep = '\t')
            gene_df['LocusPos'] = [gene_df['Pos'][i] - start_pos for i in range(len(gene_df))]
            locus_len = locus_lens[idx]
            plt.sca(axes[idx][gene_column_idx])
            axes[idx][gene_col_idx].axis('on')
            for i in range(len(gene_df)):
                gene_pos = utils.ModifyPos(gene_df['LocusPos'][i], locus_len, strands[idx])
                scale_pos = self.config.plot_scale - gene_pos / locus_len * self.config.plot_scale
            plt.plot([0, 1], [scale_pos, scale_pos], linestyle = '-', marker = 'None', color = 'black')
        plt.xlim(0, 1)
        plt.ylim(0, config.plot_scale)
        plt.xticks([], [])
        plt.yticks([], [])


class UpperTriangleUtils:
    def __init__(self, num_samples, gene_vis_utils):
        self.gene_vis_utils = gene_vis_utils
        self.num_samples = num_samples

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

class LowerTriangleUnits:
    def __init__(self, num_samples, gene_vis_utils):
        self.num_samples = num_samples
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
        return self.num_samples + self.gene_vis_utils.GetColumnWidth()

    def NumRows(self):
        return self.num_samples

    def GetLineCoordinates(self, x1, x2, y1, y2, scale):
        return [x1, x2], [scale - y1, scale - y2] 

    def SetLabels(self, axes, sample_labels):
        for idx, label in enumerate(sample_labels):
            plt.sca(axes[self.num_samples - 1, self.col_shift + idx])
            plt.xlabel(label)
