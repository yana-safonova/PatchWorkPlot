import os
import sys
import pandas as pd
import getopt

pwd = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(pwd, 'py'))

import config_utils
import visualization_utils as vis_utils
import utils
import data_utils
import tool_builder

def PrintPatchWorkLogo():
    print(' _ _ _ _ _ _')
    print('|\\  |\\ \\  |.|')
    print('|_ \\|_ _ _|_|')
    print('    |\\ \\  |.|')
    print('    |\\ \\  |.|')
    print('    |_ _ \\|_|')
    print('          |_|')

def main(command_args):
    default_params_txt = 'config.txt'
    config = config_utils.Config(default_params_txt, command_args)
    utils.PrepareDir(config.output_dir)

    input_data = data_utils.InputData(config.input_csv)
    utils.PrepareDir(config.align_dir)
    utils.PrepareDir(config.pairwise_plot_dir)

    aligner_builder = tool_builder.AlignerFactory(config)
    pairwise_aligner = aligner_builder.GetAligner()
    aligned_data = data_utils.AlignedData(input_data, pairwise_aligner, config)
    aligned_data.ReportSummaryAlignmentStats(config.align_stats_csv)
    print('Alignment stage is complete')

    print('\nVisualizing alignments...')
    visualizer_builder = tool_builder.VisualizerBuilder(config, aligned_data)
    plot_visualizer = visualizer_builder.GetPlotVisualizer()
    vis_utils.VisualizePlot(plot_visualizer, aligned_data, config)
    vis_utils.PlotPairwiseAlignments(plot_visualizer, aligned_data, config)

    print('Visualization stage is complete')

    print('\nThank you for using PatchWorkPlot!')
    PrintPatchWorkLogo()

if __name__ == '__main__':
    main(sys.argv[1:])

def cli():
    return main(sys.argv[1:])
