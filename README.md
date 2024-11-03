# PatchWorkPlot
A tool for visualization of pairwise alignments of multiple sequences as an upper or lower triangular matrix.

## Dependencies

## Input Parameters
### Mandatory parameters
`-i INPUT_CONFIG.CSV`: a configuration file containing information about input sequences in the CSV format. The configuration file contains the following colummns:
- SampleID: a unique identifier of each sequence (mandatory).
- Fasta: a complete path to each sequence in FASTA format (mandatory).
- Label: labels will be used in the output plot and, unlike SampleIDs, do not have to be unique to a sequence (mandatory).
- GeneBED: a complete path to genes in BED format (optional).  
An example of the configuration file can be found here.

`-o OUTPUR_DIR`: an output directory.

### Optional parameters
`--aligner NAME`: the name of tool used for pairwise alignment sequences. LastZ and YASS options are available. Default: `lastz`. 

`--cmap NAME`: the name of coloring map used for visualization of alignments. The minimum and maximum values of percent identity thresholds (`min-pi` and `max-pi`) will be used to determine the color of the alignment: 
- Alignments with percent identity below `min-pi` will be shown using the leftmost color in the coloring cmap.
- Alignments with percent identity above `max-pi` will be shown using the rightmost color in the coloring map.
- Alignments with percent identity between `min-pi` and `max-pi` will be projected onto the coloring map and colored accorndgly.
  
For the list of available coloring maps, please refer to [MatPlotLib documentation](https://matplotlib.org/stable/users/explain/colors/colormaps.html). Default: `Spectral`.

`--reverse-cmap BOOLEAN`: 

`--min-pi FLOAT`: 

`--max-pi FLOAT`

`--min-len INT`

`--lower`

`--show-genes`

`--transparent`

`--help`

### Default parameters


## Usage
### Visualization of IgDetective results

## Gallery
## Citation
