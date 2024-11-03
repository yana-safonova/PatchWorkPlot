# PatchWorkPlot
A tool for visualization of pairwise alignments of multiple sequences as an upper or lower triangular matrix.

## Dependencies

## Input Parameters
### Mandatory parameters
`-i INPUT_CONFIG.CSV`: a configuration file containing information about input sequences in the CSV format. The configuration file contains the following colummns:
- `SampleID`: a unique identifier of each sequence (mandatory).
- `Fasta`: a complete path to each sequence in FASTA format (mandatory).
- `Label`: labels will be used in the output plot and, unlike SampleIDs, do not have to be unique to a sequence (mandatory).
- `GeneBED`: a complete path to genes in BED format (optional).
  
An example of the configuration file can be found here.

`-o OUTPUR_DIR`: an output directory.

### Optional parameters
`--aligner NAME`: the name of tool used for pairwise alignment sequences. `lastz` (LastZ) and `yass` (YASS) options are available. Default: `lastz`. 

`--cmap NAME`: the name of coloring map used for visualization of alignments. The minimum and maximum values of percent identity thresholds (`min-pi` and `max-pi`) will be used to determine the color of the alignment: 
- Alignments with percent identity below `min-pi` will be shown using the leftmost color in the coloring cmap.
- Alignments with percent identity above `max-pi` will be shown using the rightmost color in the coloring map.
- Alignments with percent identity between `min-pi` and `max-pi` will be projected onto the coloring map and colored accorndgly.
  
For the list of available coloring maps, please refer to the [Matplotlib documentation](https://matplotlib.org/stable/users/explain/colors/colormaps.html). Default: `Spectral`.

`--reverse-cmap BOOLEAN`: if `true`, then `min-pi` and `max-pi` values of the percent identities will correspond to the rightmost and leftmost colors of the coloring map, repsectively. In case of the Spectral map, alignments with high and low percent identity will colored in red and blue, respectively. Default: `true`.  

`--min-pi FLOAT`: the alignment percent identity value that will be used to determine the color of the least similar alignments. Default: `85`.

`--max-pi FLOAT`: the alignment percent identity value that will be used to determine the color of the most similar alignments. Default: `100`.

`--min-len INT`: only alignments of lengths exceeding `min-len` will be visualized. Default: `5000`.

`--lwidth INT': the width of lines showing alignments on the final plot. Default: `1`. 

`--lower`: if specified, alignments will be visualized as a lower triangular matrix instead of an upper triangular matrix. 

`--show-genes`: if specified, gene positions will be extracted from `INPUT_CONFIG.CSV` and shown on the side of the plot. 

`--transparent`: if specified, the .PNG version of the plot will have a transparent background.  

`--help`: print help.

### Default parameters
Default values of input parameters are stored in `default_params.txt`. It covers the parameters described above as well as the parameters of LastZ alignments. The file can be modified to change the default values and avoid passing the arguments through the arguments of the command line.   

## Usage
An example of visualization of a lower-triangular patchwork plot using the `PuBuGn` coloring map (the direct orientation) with an increased line width and without gene positions: 

`python PatchWorkPlot.py -i input_config.csv -o patchwork_output --cmap PuBuGn --reverse-cmap false --lwidth 2 --lower`


Please note that the pairwise alignment is the most time-consuming step. If you want to change visualization details for previously aligned sequences, you can specify the existing output directory and change visualization parameters. E.g.: the command line:

`python PatchWorkPlot.py -i input_config.csv -o patchwork_output --cmap jet --reverse-cmap false --min-len 10000 --show-genes`

will use alignments in the `patchwork_output` direcitory and modify the patchworkplot by: 
- changing the coloring map to `jet`,
- discarding alignments shorter than 10 kbp,
- using the default parameter to draw lines, 
- showing gene positions,
- reporting an upper triangular matrix. 

### Visualization of IgDetective results

## Gallery
## Citation
