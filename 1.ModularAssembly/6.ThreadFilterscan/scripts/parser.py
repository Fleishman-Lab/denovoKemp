import pandas as pd
from scipy import special


def parse_pssm(pth: str) -> pd.DataFrame:
    df = pd.read_csv(pth, header=None, usecols=range(22), skiprows=3, delim_whitespace=True)
    with open(pth) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                break
    columns = list(line.strip().split())
    df.columns = ['POS', 'WT'] + columns[:20]
    df.set_index('POS', inplace=True)
    return df


def parse_fs(pth: str) -> pd.DataFrame:
    df = pd.read_csv(pth, delim_whitespace=True)
    df.columns = map(str.upper, df.columns)
    df.set_index('POS', inplace=True)
    return df


def parse_simple_fs(pth: str) -> pd.DataFrame:
    df = pd.read_csv(pth, delim_whitespace=True, names=['POS', 'WT', 'MUT', 'TOTAL'])
    df.set_index('POS', inplace=True)
    return df


def fs_probs(fs: pd.DataFrame, include_wt=True) -> pd.DataFrame:
    tot = fs[['WT', 'MUT', 'TOTAL']].copy()
    if include_wt:
        wt_tot = tot.groupby('POS').apply(lambda x:
                                          pd.Series([x['WT'].iloc[0], x['WT'].iloc[0], 0.0],
                                                    index=['WT', 'MUT', 'TOTAL']))
        tot = pd.concat([tot, wt_tot])
    tot['PROB'] = tot.TOTAL.groupby('POS').transform(lambda x: special.softmax(-x))
    prob = tot.pivot(columns='MUT', values='PROB').fillna(0)
    return prob
