import os
import sys
import pandas as pd
from Bio import SeqIO

sys.path.append('py')
import utils

class InputData:
    def __init__(self, data_csv):
        self.data_csv = data_csv
        self._InitiateData()

    def _InitiateData(self):
        self.data_df = pd.read_csv(self.data_csv)
        self.seq_list = []
        for i in range(len(self.data_df)):
            seq = [str(r.seq).upper() for r in SeqIO.parse(self.data_df['Fasta'][i], 'fasta')][0]
            self.seq_list.append(seq)
        self.species_names = list(self.data_df['SampleID'])
        self.locus_lens = [len(seq) for seq in self.seq_list]

    def GetLengthByIdx(self, idx):
        return self.locus_lens[idx]

    def GetSampleNameByIdx(self, idx):
        return self.data_df['SampleID'][idx]

    def GetLabelByIdx(self, idx):
        return self.data_df['Label'][idx]
    
    def GetFastaByIdx(self, idx):
        return self.data_df['Fasta'][idx]

    def NumSamples(self):
        return len(self.species_names)

    def GetGeneTableByIdx(self, idx):
        return pd.read_csv(self.data_df['GeneTxt'][idx], sep = '\t')

    #### refactor!
    def GetStartPosByIdx(self, idx):
        return self.data_df['StartPos'][idx]

class LastZPairwiseAligner:
    def __init__(self, config):
        self.config = config
        self.lastz_params = '--step=20 --notransition --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,id%'

    def AlignTwoFasta(self, fasta1, fasta2, output_fname):
        os.system('lastz ' + fasta1 + ' ' + fasta2 + ' ' + self.lastz_params + ' --output=' + output_fname)

    def GetAlignedDF(self, output_fname):
        return pd.read_csv(output_fname, sep = '\t')

class YassPairwiseAligner:
    def __init__(self, config):
        self.config = config

    def AlignTwoFasta(self, fasta1, fasta2, output_fname):
        os.system('yass -d 2 -o ' + output_fname + ' ' + fasta1 + ' ' + fasta2 + ' > /dev/null 2>&1')

    def GetAlignedDF(self, output_fname):
        lines = open(output_fname).readlines()[1:]
        df = {'id%' : [], 'start1' : [], 'end1' : [], 'start2' : [], 'end2' : []}
        for l in lines:
            splits = l.strip().split()        
            # 0name1	1name2	2id%	3alignment_length	4mismatches	5gap_opening	6start1	7end1	8start2	9end2	eval	bit_score
            df['id%'].append(splits[2] + '%')
            df['start1'].append(int(splits[6]))
            df['end1'].append(int(splits[7]))
            df['start2'].append(int(splits[8]))
            df['end2'].append(int(splits[9]))
        df = pd.DataFrame(df)
        start_list = []
        end_list = []
        strand_list = []
        for i in range(len(df)):
            if df['start2'][i] < df['end2'][i]:
                start_list.append(df['start2'][i])
                end_list.append(df['end2'][i])
                strand_list.append('+')
            else:
                start_list.append(df['end2'][i])
                end_list.append(df['start2'][i])
                strand_list.append('-')
        df['start2+'] = start_list
        df['end2+'] = end_list
        df['strand2'] = strand_list
        df['length1'] = [abs(df['start1'][i] - df['end1'][i]) for i in range(len(df))]
        df['length2'] = [abs(df['start2'][i] - df['end2'][i]) for i in range(len(df))]
        print('parsing ' + output_fname + '...')
        return df

class AlignedData:
    def __init__(self, input_data, pairwise_aligner, config):
        self.input_data = input_data
        self.pairwise_aligner = pairwise_aligner
        self.align_dir = config.align_dir
        self.config = config
        print('Computing pairwise alignments...')
        self._PerformPairwiseAlignments()
        self._ReadAlignments()
        print('Redefining strands...')
        self._RedefineStrands()
        print('Redirecting alignments...')
        self._RedirectAlignments()

    def _PerformPairwiseAlignments(self):
        self.align_dict = dict()
        for i in range(self.input_data.NumSamples()):
            print('aligning ' + self.input_data.GetSampleNameByIdx(i) + '...')

            #self dot plot
            sample_name1 = self.input_data.GetSampleNameByIdx(i)
            fasta1 = self.input_data.GetFastaByIdx(i)
            out_self = os.path.join(self.align_dir, 'self_' + str(i) + '-' + sample_name1 + '.tsv')
            if not os.path.exists(out_self):
                self.pairwise_aligner.AlignTwoFasta(fasta1, fasta1, out_self)
            print('  ' + out_self + ' was computed')
            self.align_dict[i, i] = out_self

            # pairwise dot plots
            for j in range(i + 1, self.input_data.NumSamples()):
                sample_name2 = self.input_data.GetSampleNameByIdx(j)
                fasta2 = self.input_data.GetFastaByIdx(j)
                out_pair = os.path.join(self.align_dir, 'pair_' + str(i) + '-' + sample_name1 + '_' + str(j) + '-' + sample_name2 + '.tsv')
                if not os.path.exists(out_pair):
                    self.pairwise_aligner.AlignTwoFasta(fasta1, fasta2, out_pair)
                print('  ' + out_pair + ' was computed')
                self.align_dict[i, j] = out_pair

    def _ReadAlignments(self):
        self.align_dfs = dict()
        for i in range(self.input_data.NumSamples()):
            for j in range(i, self.input_data.NumSamples()):
                raw_df = self.pairwise_aligner.GetAlignedDF(self.align_dict[i, j])
                df = raw_df.loc[(raw_df['length1'] >= self.config.min_align_len) & (raw_df['length2'] >= self.config.min_align_len)].reset_index()
                df['id%'] = [float(df['id%'][i][:-1]) for i in range(len(df))]
                self.align_dfs[i, j] = df

    def _RedefineStrands(self):
        self.strands = ['+']
        for i in range(1, self.input_data.NumSamples()):
            df = self.align_dfs[0, i] 
            max_len = max(df['length1'])
            df_max = df.loc[df['length1'] == max_len].reset_index()
            self.strands.append(df_max['strand2'][0])

    def _RedirectAlignments(self):
        for idx1, idx2 in self.align_dfs:
            len1 = self.input_data.GetLengthByIdx(idx1)
            len2 = self.input_data.GetLengthByIdx(idx2)
            df = self.align_dfs[idx1, idx2]
            directed_pos_list = []
            for i in range(len(df)):
                pos1 = df['start1'][i], df['end1'][i]
                pos2 = df['start2+'][i], df['end2+'][i]
                if df['strand2'][i] == '-':
                    pos2 = df['end2+'][i], df['start2+'][i]
                pos1 = utils.ModifyPos(pos1[0], len1, self.strands[idx1]) - 1, utils.ModifyPos(pos1[1], len1, self.strands[idx1]) - 1
                pos2 = utils.ModifyPos(pos2[0], len2, self.strands[idx2]) - 1, utils.ModifyPos(pos2[1], len2, self.strands[idx2]) - 1
                directed_pos_list.append([pos1[0], pos1[1], pos2[0], pos2[1]])
            df['start1_dir'] = [p[0] for p in directed_pos_list]
            df['end1_dir'] = [p[1] for p in directed_pos_list]
            df['start2_dir'] = [p[2] for p in directed_pos_list]
            df['end2_dir'] = [p[3] for p in directed_pos_list]

    def GetLengthByIdx(self, idx):
        return self.input_data.GetLengthByIdx(idx)

    def GetSampleNameByIdx(self, idx):
        return self.input_data.GetSampleNameByIdx(idx)

    def GetLabelByIdx(self, idx):
        return self.input_data.GetLabelByIdx(idx)

    def GetFastaByIdx(self, idx):
        return self.input_data.GetFastaByIdx(idx)

    def NumSamples(self):
        return self.input_data.NumSamples()

    def GetGeneTableByIdx(self, idx):
        return self.input_data.GetGeneTableByIdx(idx)

    #### refactor!
    def GetStartPosByIdx(self, idx):
        return self.input_data.GetStartPosByIdx(idx)

    def GetAlignmentDF(self, idx1, idx2):
        return self.align_dfs[idx1, idx2]

    def GetStrandByIdx(self, idx):
        return self.strands[idx]

    def ReportSummaryAlignmentStats(self, output_fname):
        stats_df = {'Label1' : [], 'Label2' : [], 'Idx1' : [], 'Idx2' : [], 'PI' : []}
        for idx1, idx2 in self.align_dfs:
            if idx1 == idx2:
                continue
            df = self.align_dfs[idx1, idx2]
            for i in range(len(df)):
                stats_df['Label1'].append(self.GetSampleNameByIdx(idx1))
                stats_df['Idx1'].append(idx1)
                stats_df['Label2'].append(self.GetSampleNameByIdx(idx2))
                stats_df['Idx2'].append(idx2)
                stats_df['PI'].append(df['id%'][i])
        stats_df = pd.DataFrame(stats_df)
        stats_df.to_csv(output_fname, index = False)

    def IndexPairIterator(self):
        for idx1, idx2 in self.align_dfs:
            yield idx1, idx2
