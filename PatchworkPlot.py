import os
import sys
import pandas as pd
from Bio import SeqIO

sys.path.append('py')
import visualization_utils as vis_utils
import utils
import data_utils

class Config:
    def __init__(self, output_dir):
        #### alignment params
        self.min_align_len = 10000
        self.alignment_method = 'lastz' # or 'yass'  

        #### visualization params
        self.pi_min = 90
        self.pi_max = 100
        self.plot_scale = 1000
        self.cmap = 'Spectral'
        self.cmap_reverse = True
        self.upper_triangle = True
        self.linewidth = 1
        self.show_genes = True

        #### output params
        self.transparent = False
        self.output_dir = output_dir
        self.align_dir = os.path.join(self.output_dir, 'pairwise_alignments')

def main(data_csv, output_dir):
    config = Config(output_dir)
    utils.PrepareDir(output_dir)

    input_data = data_utils.InputData(data_csv)
    utils.PrepareDir(config.align_dir)
    pairwise_aligner = data_utils.LastZPairwiseAligner(config.align_dir)
    aligned_data = data_utils.AlignedData(input_data, pairwise_aligner, config)
    aligned_data.ReportSummaryAlignmentStats(os.path.join(output_dir, 'pi_stats.csv')) 
    print('Alignment stage is complete')

    print('\nVisualizing alignments...')
    gene_visualizer = vis_utils.SimpleGeneVisualizer(config, aligned_data)
    plot_visualizer = vis_utils.LowerTriangleUtils(aligned_data, gene_visualizer)
    vis_utils.VisualizePlot(plot_visualizer, aligned_data, config)
    print('Visualization stage is complete')

    print('\nThank you for using PatchworkPlot!')

if __name__ == '__main__':
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    main(input_csv, output_dir)
