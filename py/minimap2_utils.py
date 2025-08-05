import pandas as pd

class Minimap2Reader:

    def __init__(self, path):
        self.path = path

    def ParsePaf(self):
        df = pd.read_csv(self.path, header=None, sep='\t',
                         usecols=range(12), comment="#", engine="python")

        df = df.rename(columns={5: '#name1',
                                4: 'strand1',
                                7: 'start1',
                                8: 'end1',
                                0: 'name2',
                                2: 'start2+',
                                3: 'end2+',
                                9: 'id%'})
        df['strand2'] = '+'

        df['length1'] = (df['end1'].astype(int) - df['start1'].astype(int))
        df['length2'] = (df['end2+'].astype(int) - df['start2+'].astype(int))

        df['id%'] = df['id%'] / df[['length1', 'length2']].min(axis=1) * 100
        df['id%'] = df['id%'].apply(lambda x: str(round(x, 1)) + '%')

        df = df[['#name1', 'strand1', 'start1', 'end1', 'length1', 'name2',
                 'strand2', 'start2+', 'end2+', 'length2', 'id%']]

        return df
