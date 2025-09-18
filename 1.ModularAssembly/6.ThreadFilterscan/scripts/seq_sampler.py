from typing import List, Tuple, Union

import numpy as np
import pandas as pd


def sample_indices(adj: np.array, indices: List[int]) -> List[int]:
    sampled_indices: List[int] = []
    while len(indices):
        sam_idx = np.random.choice(indices)
        sampled_indices.append(sam_idx)
        neigh = np.flatnonzero(adj[sam_idx])
        indices = list(set(indices).difference(neigh))

    return sampled_indices


def sample_seq(adj: np.array, probs: pd.DataFrame) -> List[Tuple[int, str]]:
    """
    The assumption here is that position is sequential and starts at 1, that is the
    first AA is numbered as 1.
    The random sample goes as follows:
        As long there are positions we can sample repeat:
            1. Randomly select a position
            2. Remove adjacent positions
        Now that we have positions - for each position sample a letter from probs.
    :param adj:
    :param probs:
    :return:
    """

    sampled_indices = np.array(sample_indices(adj, (probs.index - 1).values.tolist())) + 1
    return sorted([(idx, np.random.choice(probs.columns, p=probs.loc[idx])) for idx in sampled_indices])


def greedy_seq(adj: np.array, fs_table: pd.DataFrame) -> List[Tuple[int, str]]:
    """
    The assumption here is that position is sequential and starts at 1, that is the
    first AA is numbered as 1.
    The random sample goes as follows:
        As long there are positions we can sample repeat:
            1. Randomly select a position
            2. Remove adjacent positions
        Now that we have positions - for each position sample a letter from probs.
    :param adj:
    :param probs:
    :return:
    """
    lowest_energy = fs_table[fs_table.TOTAL < 0].groupby('POS').apply(
        lambda x: x[x['TOTAL'] == x['TOTAL'].min()]).reset_index(1, drop=True)
    indices = lowest_energy.index.values
    select_indices: List[Tuple[int, str]] = []
    while len(indices):
        index = lowest_energy.loc[indices, 'TOTAL'].idxmin()
        mut = lowest_energy.loc[index, 'MUT']
        select_indices.append((index, mut))
        neigh = np.flatnonzero(adj[index - 1]) + 1
        indices = np.array(list(set(indices).difference(neigh)))
    return sorted(select_indices)


def _seqlist(seq: Union[pd.DataFrame, str, List[str]]) -> List[str]:
    if isinstance(seq, pd.DataFrame):
        return seq['WT'].values.tolist()
    elif isinstance(seq, str):
        return list(seq)
    elif isinstance(seq, list):
        return seq
    else:
        raise ValueError("Provided sequence is not pandas/list/str")


def get_mutations(sample: List[Tuple[int, str]], wt_seq: Union[pd.DataFrame, str, List[str]]) -> List[Tuple[int, str]]:
    """
    wt_seq can be either a string, list of AAs or a PSSM dataframe
    :param sample:
    :param wt_seq:
    :return:
    """
    seq = _seqlist(wt_seq)
    return [(idx, aa) for (idx, aa) in sample if aa != seq[idx - 1]]


def mutated_seq(mutations: List[Tuple[int, str]], wt_seq: Union[pd.DataFrame, str, List[str]]) -> str:
    seq = _seqlist(wt_seq)
    for idx, aa in mutations:
        seq[idx - 1] = aa
    return ''.join(seq)


def sequence_score(seq: Union[pd.DataFrame, str, List[str]], wt_seq: Union[pd.DataFrame, str, List[str]],
                   fs_table: pd.DataFrame) -> List[int]:
    seq = _seqlist(seq)
    wt_seq = _seqlist(wt_seq)
    scores = []
    for i, (seq_aa, wt_aa) in enumerate(zip(seq, wt_seq)):
        if seq_aa == wt_aa:
            scores.append(0.0)
            continue
        scores.append(fs_table[(fs_table.index == i + 1) & (fs_table.WT == wt_aa.upper()) & (
                fs_table.MUT == seq_aa.upper())].TOTAL.item())
    return scores
