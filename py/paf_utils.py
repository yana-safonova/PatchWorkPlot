import pandas as pd

class PafReader:

    def __init__(self, path, sep):
        self.path = path
        self.sep = sep

    def ParsePaf(self):

        df = pd.read_csv(self.path, header=None, sep=self.sep,
                         usecols=range(10), comment="#", engine="python")

        df = df.rename(columns={5: '#name1',
                                4: 'strand1',
                                7: 'start1',
                                8: 'end1',
                                0: 'name2',
                                2: 'start2+',
                                3: 'end2+',
                                9: 'id%'})

#        df['strand2'] = '+'

        df['length1'] = (df['end1'].astype(int) - df['start1'].astype(int))
        df['length2'] = (df['end2+'].astype(int) - df['start2+'].astype(int))

#        mask = df['strand1'] == '-'
#        df.loc[mask, ['start2+', 'end2+']] = df.loc[mask, ['end2+', 'start2+']].values
        df['strand2'] = df['strand1']
        df['strand1'] = '+'

        if self.sep == '\t':
            df['id%'] = df['id%'] / df[['length1', 'length2']].min(axis=1) * 100
        df['id%'] = df['id%'].apply(lambda x: str(round(x, 1)) + '%')

        df = df[['#name1', 'strand1', 'start1', 'end1', 'length1', 'name2',
                 'strand2', 'start2+', 'end2+', 'length2', 'id%']]

        return df
