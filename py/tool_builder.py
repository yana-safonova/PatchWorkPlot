import os
import sys

import data_utils
import visualization_utils as vis_utils

class AlignerFactory:
    def __init__(self, config):
        self.config = config

    def GetAligner(self):
        if self.config.alignment_method == 'yass':
            return data_utils.YassPairwiseAligner(self.config)
        elif self.config.alignment_method == 'minimap2':
            return data_utils.Minimap2Aligner(self.config)
        elif self.config.alignment_method == 'mashmap':
            return data_utils.MashmapAligner(self.config)
        elif self.config.alignment_method == 'custom':
            return data_utils.CustomAligner(self.config)
        return data_utils.LastZPairwiseAligner(self.config)


class VisualizerBuilder:
    def __init__(self, config, aligned_data):
        self.config = config
        self.aligned_data = aligned_data

    def _GetGeneVisualizer(self):
        if self.config.show_annotation:
            return vis_utils.SimpleGeneVisualizer(self.config, self.aligned_data)
        return vis_utils.EmptyGeneVisualizer(self.config, self.aligned_data)

    def GetPlotVisualizer(self):
        gene_visualizer = self._GetGeneVisualizer()
        if self.config.upper_triangle:
            return vis_utils.UpperTriangleUtils(self.aligned_data, gene_visualizer)
        return vis_utils.LowerTriangleUtils(self.aligned_data, gene_visualizer)
        
