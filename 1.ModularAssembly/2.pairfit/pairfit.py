#!/usr/bin/env python3

import os
import logging
import argparse
import sys
from flab.utils.parse_actions import FullPaths, FullPathsList


__author__ = 'Rosalie Lipsh-Sokolik'


LGR = logging.getLogger(__file__)


def get_pdb_name(pdb_path):
    return os.path.basename(pdb_path).replace('.pdb', '').replace('.gz', '')


def load_pdb(path, name=None):
    if not name:
        name = get_pdb_name(path)
    cmd.load(path, object=name)
    cmd.remove("het")


def get_resi_list(pdb):
    space = {'residues': list()}
    cmd.iterate(pdb, 'residues.append(resi)', space=space)
    residues = sorted([int(r) for r in set(space['residues'])])
    return residues


def make_stems_selection(stem1_index, stem2_index, pdb, strech=3):
    """
    :param stem1_index: the very start of the fragment
    :param stem2_index: the last residue of the fragment
    :param strech: length of the stem
    """
    selection = '{} and resi {}-{} and name c+ca+n+o'
    name1 = '{}_stem1'.format(pdb)
    cmd.select(name1, selection.format(pdb, stem1_index, stem1_index+strech-1))
    name2 = '{}_stem2'.format(pdb)
    cmd.select(name2, selection.format(pdb, stem2_index-strech+1, stem2_index))
    return (name1, name2)


def is_continues(residues):
    """Makes sure a strech of 3 resiues is continues
    :param residues: ordered list of resi
    :return: True is the resi are continues
    """
    for i in range(len(residues)-1):
        if residues[i+1] - residues[i] != 1:
            return False
    return True    


def pair_fit(target_stem1, target_stem2, template_stem1, template_stem2):
    """All params are names of selections. The target will move on the PyMOL
    session
    """
    score = cmd.pair_fit('{} + {}'.format(target_stem1, target_stem2),
                         '{} + {}'.format(template_stem1, template_stem2))
    return score


def save_structure(pdb, path=None):
    """
    :param pdb: the name of the pymol object to save
    :param path: where to save it
    """
    if not path:
        path = '{}_aligned.pdb'.format(get_pdb_name(pdb))
    cmd.save(path, pdb)
    return path



# Main
strech = 3

data  = {'template': dict(), 'target': dict()}
data['template']['path'] = sys.argv[-2]
data['target']['path'] = sys.argv[-1]

for fragment, d in data.items():
    load_pdb(d['path'], fragment)
    
    resi = get_resi_list(fragment)
    d['stem1'] = resi[0]
    if not is_continues(resi[:strech]):
        msg = 'Start stem of {} is not contenues'.format(d["path"])
        raise ValueError(msg)
    d['stem2'] = resi[-1]
    if not is_continues(resi[(-1*strech):]):
        msg = 'End stem of {}is not contenues'.format(d["path"])
        raise ValueError(msg)

    d['stem1_sel'], d['stem2_sel'] = make_stems_selection(d['stem1'],
                                                          d['stem2'], 
                                                          fragment)
pair_fit(data['target']['stem1_sel'], data['target']['stem2_sel'],
         data['template']['stem1_sel'], data['template']['stem2_sel'])
index = data['target']['path'].rfind('.pdb')
save_structure('target',
               path='{}_aligned.pdb'.format(data['target']['path'][:index]))
