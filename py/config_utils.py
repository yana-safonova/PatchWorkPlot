import os
import sys
import getopt

class Config:
    def __init__(self, default_txt, command_args):
        self._ReadDefaultParameters(default_txt)
        self._ParseCommandLineParams(command_args)
        self._CheckParameterCompleteness()

    def _ReadDefaultParameters(self, default_txt):
        #### input csv
        self.input_csv = ''

        #### alignment params
        self.min_align_len = 5000
        self.alignment_method = 'lastz' # or 'yass', 'minimap2', 'mashmap', 'custom'
        self.lastz_params = '--step=20 --notransition'
        self.minimap2_params = '--secondary=yes -P -k 10 -w 5 --no-long-join -r 100 -g 50'
        self.mashmap_params = '--pi 70'

        #### visualization params
        self.pi_min = 85
        self.pi_max = 100
        self.plot_scale = 1000
        self.cmap = 'Spectral'
        self.cmap_reverse = True
        self.color = ''
        self.upper_triangle = True
        self.linewidth = 1
        self.show_annotation = False
        self.hide_legend = False
        self.show_breakpoints = False
        self.bp_color = '#7F7F7F'
        self.bp_min_len = 10000
        self.bp_linewidth = 0.2

        #### output params
        self.transparent = False
        self.output_dir = ''
        self.verbose = 2

    def _ParseCommandLineParams(self, command_args):
        opts = []
        try:
            opts, args = getopt.getopt(command_args, 'i:o:hv::',  ['min-pi=', 'max-pi=', 'aligner=',
                                                               'minimap2-params=', 'mashmap-params='
                                                               'min-len=', 'cmap=', 'reverse-cmap=',
                                                               'color=', 'lower', 'lwidth=', 'show-annot', 
                                                               'show-bp', 'bp-color=', 'bp-min-len=', 'bp-lwidth=',
                                                               'transparent', 'help',
                                                               'verbose=', 'v=', 'hide-legend'])
        except:
            print('Error')
        for opt, arg in opts:
            if opt == '-i':
                self.input_csv = arg
            elif opt == '-o':
                self.output_dir = arg
            elif opt == '--min-pi':
                self.pi_min = float(arg)
            elif opt == '--max-pi':
                self.pi_max =  float(arg)
            elif opt == '--aligner':
                self.alignment_method = arg
            elif opt == '--minimap2-params':
                self.minimap2_params = arg
            elif opt == '--mashmap-params':
                self.mashmap_params = arg
            elif opt == '--min-len':
                self.min_align_len = int(arg)
            elif opt == '--cmap':
                self.cmap = arg
            elif opt == '--reverse-cmap':
                self.cmap_reverse = arg.lower() == 'true'
            elif opt == '--color':
                self.color = arg
            elif opt == '--lower':
                self.upper_triangle = False
            elif opt == '--lwidth':
                self.linewidth = float(arg)
            elif opt == '--show-annot':
                self.show_annotation = True
            elif opt == '--show-bp':
                self.show_breakpoints = True
            elif opt == '--bp-color':
                self.bp_color = arg
            elif opt == '--bp-min-len':
                self.bp_min_len = int(arg)
            elif opt == '--bp-lwidth':
                self.bp_linewidth = float(arg)
            elif opt == '--transparent':
                self.transparent = True
            elif opt == '--hide-legend':
                self.hide_legend = True
            elif opt == '--verbose' or opt == '-v':
                self.verbose = int(arg)
            elif opt == '--help' or opt == '-h':
                self.PrintHelpMessage()
                sys.exit(0)

    def _CheckParameterCompleteness(self):
        if self.input_csv == '' or not os.path.exists(self.input_csv):
            print('ERROR: config file \"' + self.input_csv + '\" (-i) was empty or not found')
            sys.exit(1)
        if self.output_dir == '':
            print('ERROR: output directory (-o) was not speficied')
            sys.exit(1)

        self.align_dir = os.path.join(self.output_dir, 'pairwise_alignments')
        self.align_stats_csv = os.path.join(self.output_dir, 'alignment_stats.csv')
        self.pairwise_plot_dir = os.path.join(self.output_dir, 'pairwise_dotplots')

    def PrintHelpMessage(self):
        print('python PatchWorkPlot.py -i INPUT_CONFIG.CSV -o OUTPUT_DIR [--aligner METHOD --min-pi FLOAT --max-pi FLOAT --min-len INT --cmap CMAP_NAME --reverse-cmap BOOL --lower --lwidth INT --show-annot --transparent --hide-legend --verbose/-v VERBOSE --help/-h]')
