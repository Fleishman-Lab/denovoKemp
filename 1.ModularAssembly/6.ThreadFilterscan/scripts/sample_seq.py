import fire
import numpy as np
from itertools import islice
from scipy.spatial import distance

from pdb_adj import res_dist
from parser import fs_probs, parse_fs, parse_pssm, parse_simple_fs
from seq_sampler import mutated_seq, sample_seq, greedy_seq


def sample(pssm: str, pdb: str, fs_table: str, dist_cutoff: int = 8, include_wt=False):
    """
    Sample a sequence from non-overlapping positions with fs scores as probabilities.
    :param pssm: Path to PSSM file
    :param pdb: Path to PDB file
    :param fs_table: Path to FS results table
    :param dist_cutoff: Under which distance do positions are considered adjacent
    :param include_wt: Whether to indlude WT in the sample probabilities which energy 0.0
    :return: None
    """
    pssm = parse_pssm(pssm)
    dist_mat = distance.squareform(res_dist(pdb))
    adj = (dist_mat < dist_cutoff).astype(np.int)
    probs = fs_probs(parse_simple_fs(fs_table), include_wt=include_wt)
    seq = sample_seq(adj, probs)
    print(seq)
    print(mutated_seq(seq, pssm))


def greedy(pssm: str, pdb: str, fs_table: str, dist_cutoff: int = 8):
    """
    Greedily select mutations by FS energy scores.
    After a mutation is selected adjacent positions are not allowed to change.
    :param pssm: Path to PSSM file or sequence string
    :param pdb: Path to PDB file
    :param fs_table: Path to FS results table
    :param dist_cutoff: Under which distance do positions are considered adjacent
    :return: None
    """
    try:
        pssm = parse_pssm(pssm)
    except OSError:
        pass

    dist_mat = distance.squareform(res_dist(pdb))
    adj = (dist_mat < dist_cutoff).astype(np.int)
    fs_table = parse_simple_fs(fs_table)
    seq = greedy_seq(adj, fs_table)
    print(seq)
    print(mutated_seq(seq, pssm))


def mut_range(wt: str, mut: str, chain: str = 'A', start: int = None, stop: int = None) -> str:
    """
    Create a mut-range to be supplied to filterscan_create_run.
    Mutations are held fixed while unmutated positions are allowed to be scanned.
    :param wt: str The WT sequence
    :param mut: str The mutant sequence
    :param chain: Which chain to add to mut-range
    :param stop: At which position to start - uses python numbering - i.e first is 0
    :param start: At which position to end - uses python numbering - i.e first is 0
    :return: str mut-range string
    """
    chain = chain.upper()
    mut_range = []
    open_range = False
    for i, (wt_aa, mut_aa) in islice(enumerate(zip(wt.upper(), mut.upper())), start, stop):
        """
        We need to take care of 4 cases:
        1. A range is not open and wt_aa == mut_aa - Need to open a new range
        2. A range is not open and wt_aa != mut_aa - Skip the current letter
        3. A range is open and wt_aa == mut_aa - Extend the current range (skip)
        4. A range is open and wt_aa != mut_aa - Close the current range
        """
        if open_range:
            if wt_aa == mut_aa:
                continue
            else:
                mut_range.append(f'{i + 1}{chain},')
                open_range = False
        else:
            if wt_aa == mut_aa:
                mut_range.append(f'{i + 1}{chain}:')
                open_range = True
            else:
                continue
    if open_range:
        mut_range.append(f'{i + 2}{chain}')
    return ''.join(mut_range)


if __name__ == '__main__':
    # See: https://github.com/google/python-fire/issues/188#issuecomment-631419585
    def Display(lines, out):
        text = "\n".join(lines) + "\n"
        out.write(text)


    from fire import core

    core.Display = Display

    fire.Fire({
        "sample": sample,
        "greedy": greedy,
        "mut-range": mut_range,
    })
