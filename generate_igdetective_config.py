import os
import sys
import pandas as pd

pwd = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(pwd, 'py'))

import utils

igdetect_dirs = sys.argv[1].split()
locus = sys.argv[2]
output_dir = sys.argv[3]

utils.PrepareDir(output_dir)

df = {'SampleID' : [], 'Label' : [], 'Fasta' : [], 'Annotation' : [], 'Strand' : []}
for igdetect_dir in igdetect_dirs:
    locus_dir = os.path.join(igdetect_dir, 'refined_ig_loci')
    summary_csv = os.path.join(locus_dir, 'summary.csv')
    locus_seq_dir = os.path.join(locus_dir, 'igloci_fasta')
    summary_df = pd.read_csv(summary_csv)
    locus_df = summary_df.loc[summary_df['Locus'] == locus].reset_index()
    if len(locus_df) != 1:
        continue

    #### extracting FASTA and gene tables
    label = '_'.join(os.path.basename(igdetect_dir).split('_')[:-1])
    fasta_fname = locus + '_' + locus_df['Contig'][0] + '_' + str(locus_df['NumV'][0]) + 'Vs.fasta'
    fasta_fname = os.path.abspath(os.path.join(locus_seq_dir, fasta_fname)) 
    igdetect_txt = os.path.join(igdetect_dir, 'combined_genes_' + locus + '.txt')
    igdetect_df = pd.read_csv(igdetect_txt, sep = '\t')

    #### converting the gene table to BED
    annot_df = {'Start' : [], 'End' : [], 'Color' : []}
    locus_start = locus_df['StartPos'][0]
    for i in range(len(igdetect_df)):
        gene_start_pos = igdetect_df['Pos'][i] - locus_start
        annot_df['Start'].append(gene_start_pos)
        annot_df['End'].append(gene_start_pos + len(igdetect_df['Sequence'][i]))
        annot_df['Color'].append('0,0,0')
    annot_df = pd.DataFrame(annot_df)
    annot_bed = os.path.abspath(os.path.join(output_dir, label + '.bed'))
    annot_fh = open(annot_bed, 'w')
    for i in range(len(annot_df)):
        annot_fh.write('NA\t' + str(annot_df['Start'][i]) + '\t' + str(annot_df['End'][i]) + '\t' + '\t'.join(['NA'] * 5) + '\t' + annot_df['Color'][i] + '\n')
    annot_fh.close()

    #### updating config DF
    df['Label'].append(label)
    df['SampleID'].append(label)
    df['Fasta'].append(fasta_fname)
    df['Annotation'].append(annot_bed)
    df['Strand'].append('')

df = pd.DataFrame(df)
output_csv = os.path.join(output_dir, 'config.csv')
df.to_csv(output_csv, index = False)
