import os
import sys
import pandas as pd
import py.paf_utils

pwd = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(pwd, 'py'))

#### python3 generate_alignment_files.py {alignment.config.tsv} {output} || name 0-x, 1-y

align_config_path = sys.argv[1].split()[0]
out_dir = sys.argv[2].split()[0]

sep = '\t' if align_config_path.endswith('.tsv') else ','
align_config = pd.read_csv(align_config_path, sep=sep)

for index, row in align_config.iterrows():

    name_1 = row['name1']
    name_2 = row['name2']
    paf_path = row['pafPath']

    print(f'Generating alignment file from {paf_path}...')
    parser = py.paf_utils.PafReader(paf_path, sep='\t')
    paf = parser.ParsePaf()

    if name_1 == name_2:
        paf.to_csv(f'{out_dir}/self_{name_1}.tsv', sep='\t', index=False)
    else:
        paf.to_csv(f'{out_dir}/pair_{name_1}_{name_2}.tsv', sep='\t', index=False)