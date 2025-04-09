import functools
from itertools import combinations
from typing import Callable, Iterator

import numpy as np
from Bio import PDB
from Bio.PDB import Atom, Vector
from scipy.spatial import distance

parser = PDB.PDBParser(QUIET=True)

backbone_atoms = [' N  ', ' CA ', ' C  ', ' O  ']
gly_backbone_atoms = [' N  ', ' C  ', ' O  ']


def filter_hydrogen(atoms: Iterator[Atom.Atom]) -> Iterator[Atom.Atom]:
    return filter(lambda x: 'H' not in x.get_fullname(), atoms)


def filter_backbone(atoms: Iterator[Atom.Atom]) -> Iterator[Atom.Atom]:
    return filter(lambda x: x.get_fullname() not in backbone_atoms, atoms)


@functools.lru_cache(maxsize=None)
def res_atom_mat(atoms: Iterator[Atom.Atom]) -> np.array:
    return np.vstack(list(map(Vector.get_array, map(Atom.Atom.get_vector, atoms))))


def no_bb_dist(op: Callable) -> Callable:
    def dist_m(res1_atoms: Iterator[Atom.Atom], res2_atoms: Iterator[Atom.Atom]) -> float:
        return op(
            [
                op(distance.cdist(res_atom_mat(filter_backbone(res1_atoms)), res_atom_mat(res2_atoms))),
                op(distance.cdist(res_atom_mat(filter_backbone(res2_atoms)), res_atom_mat(res1_atoms)))
            ]
        )

    return dist_m


def w_bb_dist(op: Callable) -> Callable:
    def dist_m(res1_atoms: Iterator[Atom.Atom], res2_atoms: Iterator[Atom.Atom]) -> float:
        return op(distance.cdist(res_atom_mat(res1_atoms), res_atom_mat(res2_atoms)))

    return dist_m


def no_h_dist(dist: Callable) -> Callable:
    def no_h_dist_m(res1_atoms: Iterator[Atom.Atom], res2_atoms: Iterator[Atom.Atom]) -> float:
        return dist(filter_hydrogen(res1_atoms), filter_hydrogen(res2_atoms))

    return no_h_dist_m


def res_dist(pdb: str, op: Callable = np.min, filter_bb: bool = False, filter_h: bool = False) -> np.array:
    struct = parser.get_structure('A', pdb)
    res_iter1 = struct.get_residues()
    dist_op: Callable = w_bb_dist(op)
    if filter_bb:
        dist_op = no_bb_dist(op)
    if filter_h:
        dist_op = no_h_dist(dist_op)
    dist = [dist_op(res1, res2) for res1, res2 in combinations(res_iter1, 2)]
    return dist
