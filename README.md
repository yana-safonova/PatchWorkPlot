# PatchWorkPlot
A tool for visualization of pairwise alignments of multiple sequences as [dot plots](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)) arranged into an upper or lower triangular matrix.

## Dependencies
- [Pandas](https://anaconda.org/anaconda/pandas)
- [BioPython](https://anaconda.org/conda-forge/biopython)
- [LASTZ](https://anaconda.org/bioconda/lastz) and/or [YASS](https://anaconda.org/bioconda/yass)

## Input Parameters
### Required parameters
`-i INPUT_CONFIG.CSV`: a configuration file containing information about input sequences in the CSV format. The configuration file contains the following colummns:
- `SampleID`: a unique identifier of each sequence (required).
- `Fasta`: a complete path to each sequence in FASTA format (required).
- `Label`: labels will be used in the output plot and, unlike SampleIDs, do not have to be unique to a sequence and can be empty (required).
- `Annotation`: a complete path to annotation in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (optional).
  
An example of the configuration file can be found here.

`-o OUTPUR_DIR`: the name of output directory. If the directory does not exist, it will be created.

### Optional parameters
`--aligner NAME`: the name of tool used for pairwise alignment sequences. `lastz` (LASTZ) and `yass` (YASS) options are available. Default: `lastz`. 

`--cmap NAME`: the name of coloring map used for visualization of alignments. The minimum and maximum values of percent identity thresholds (`min-pi` and `max-pi`) will be used to determine the color of the alignment: 
- Alignments with percent identity below `min-pi` will be shown using the leftmost color in the coloring cmap.
- Alignments with percent identity above `max-pi` will be shown using the rightmost color in the coloring map.
- Alignments with percent identity between `min-pi` and `max-pi` will be projected onto the coloring map and colored accorndgly.
  
For the list of available coloring maps, please refer to the [Matplotlib documentation](https://matplotlib.org/stable/users/explain/colors/colormaps.html). Default: `Spectral`.

`--reverse-cmap BOOLEAN`: if `true`, then `min-pi` and `max-pi` values of the percent identities will correspond to the rightmost and leftmost colors of the coloring map, respectively. In case of the Spectral map, alignments with high and low percent identity will colored in red and blue, respectively. Default: `true`.  

`--min-pi FLOAT`: the alignment percent identity value that will be used to determine the color of the least similar alignments. Default: `85`.

`--max-pi FLOAT`: the alignment percent identity value that will be used to determine the color of the most similar alignments. Default: `100`.

`--min-len INT`: only alignments of lengths exceeding `min-len` will be visualized. Default: `5000`.

`--lwidth INT`: the width of lines showing alignments on the final plot. Default: `1`. 

`--lower`: if specified, alignments will be visualized as a lower triangular matrix instead of an upper triangular matrix. 

`--show-annot`: if specified, annotations will be extracted from `INPUT_CONFIG.CSV` (column `Annotation`) and shown on the side of the plot. 

`--transparent`: if specified, the .PNG version of the plot will have a transparent background.  

`--help`: print help.

### Default values of parameters
Default values of input parameters are stored in `default_params.txt`. It covers the parameters described above as well as the parameters of LASTZ alignments. The file can be modified to change the default values and avoid passing the arguments through the arguments of the command line.   

## Usage
An example of visualization of a lower-triangular patchwork plot using the `PuBuGn` coloring map (the direct orientation) with an increased line width and without gene positions: 

`python PatchWorkPlot.py -i input_config.csv -o patchwork_output --cmap PuBuGn --reverse-cmap false --lwidth 2 --lower`


Please note that the pairwise alignment is the most time-consuming step. If you want to change visualization of previously aligned sequences, you can rerun PatchWorkPlot specifying the existing output directory through `-o` and changing the desired visualization parameters. E.g.: the command line:

`python PatchWorkPlot.py -i input_config.csv -o patchwork_output --cmap jet --reverse-cmap false --min-len 10000 --show-annot`

will use alignments in the `patchwork_output` direcitory and modify the patchworkplot by: 
- changing the coloring map to `jet`,
- discarding alignments shorter than 10 kbp,
- using the default parameter to draw lines, 
- showing gene positions,
- reporting an upper triangular matrix. 

### Visualization of IgDetective results
PatchWorkPlot is useful to visualize highly repetitive sequences or sequences with a high density of structural variations such as immunoglobulin (IG) and T-cell receptor (TCR) loci. A script `generate_igdetective_config.py` simplifies analysis of adaptive immune loci annotated using the IgDetective tool and generates a config file that can be used as an input to PatchWorkPlot. To run the script, use the following command line:

`python generate_igdetective_config.py PATHS_TO_IGDETECTIVE_DIRS LOCUS OUTPUT_DIR`

where:
- `PATHS_TO_IGDETECTIVE_DIRS` is a space- or comma-separated paths to output directories of IgDetective. If the paths are separated by spaces, make sure to put them in double quotes: `"PATH_1 PATH_2 ... PATH_N"`.
- `LOCUS`: a type of adaptive immune locus for which the config will be generated. Available options are `IGH, IGK, IGL, TRA, TRB, TRG`.
- `OUTPUT_DIR`: the name of output directory. If the directory does not exist, it will be created.

#### Example of joint usage of IgDetective & PatchWorkPlot 
The directory [`test_dataset`](test_dataset) includes five IgDetective directories containing results of IG/TR locus annotation for five cat genomes: 
- `test_dataset/01mPumCon1.1_hap1_igdetective`: the mountain lion (_Puma concolor_), haplotype 1, accession: 
- `test_dataset/02mPumCon1.1_hap2_igdetective`: the mountain lion (_Puma concolor_), haplotype 2, accession: 
- `test_dataset/03mNeoNeb1_igdetective`: the clouded leopard (_Neofelis nebulosa_), accession: 
- `test_dataset/04mLynRuf1_igdetective`: the bobcat (_Lynx rufus_), accession: .
- `test_dataset/05mFelCat1_igdetective`: the domestic cat (_Felis catus_), accession: .

The following command lines generates a configuration file for immunoglobulin heavy chain (IGH) loci and converts IGH gene files to BED format:

`python generate_igdetective_config.py "test_dataset/01mPumCon1.1_hap1_igdetective test_dataset/02mPumCon1.1_hap2_igdetective test_dataset/03mNeoNeb1_igdetective test_dataset/04mLynRuf1_igdetective test_dataset/05mFelCat1_igdetective" IGH cats_IGH_configuration`

Then, PatchWorkPlot takes the compiled configuration file and visualizes pairwise alignments of the IGH loci. The `--show-annot` option is used to illustrate positions of IGH genes predicted by IgDetective: 

`python PatchWorkPlot.py -o cats_IGH_configuration/config.csv -o cats_IGH_patchworkplot --show-annot`

## Gallery
| Annotation | Upper triangle | Lower triangle |
| ----| ------| ------|
| No annotation | <img src="examples/upper_no_annot.png" alt="upper_no_annotation" width="400"/> | <img src="examples/lower_no_annot.png" alt="lower_no_annotation" width="400"/>|
| Parameters | `default` | `--lwidth 2 --min-len 20000 --lower --min-pi 90 --max-pi 95`|
| With annotation | <img src="examples/upper_annot.png" alt="upper_annotation" width="400"/> | <img src="examples/lower_annot.png" alt="lower_annotation" width="400"/> |
| Parameters | `--show-annot --lwidth 2 --min-len 10000 --cmap Greens --reverse-cmap false --min-pi 80` | `--show-annot --lower --lwidth 2 --min-len 10000 --cmap viridis --min-pi 75` |

## Citation
To be added
