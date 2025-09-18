"""
scripts and classes to use PyMOL
"""

import __main__

__main__.pymol_argv = ["pymol", "-qc"]  # Pymol: quiet and no GUI
import pymol

pymol.finish_launching()

import os
from typing import Dict, List
from pymol import cmd
from collections import OrderedDict


def load_protein(prot):
    """
    Loads prot to Pymol.
    :param prot:
    :return: the name of the protein, -1 if the file doesn't exists
    """
    loaded = cmd.get_object_list()
    if prot in loaded:
        return prot
    elif os.path.isfile(prot):
        name = os.path.basename(prot).replace(".pdb", "")
        cmd.load(prot)
        return name
    else:
        return -1


def order_alignment(alignment, first):
    """
    Orders alignment object to have the protein named "first" to be in index 0
    of every item
    :param alignment: list of lists of tuples  [(prot_a, index), (prot_b, index)]
    :param first: name of protein to be in index 0
    :return: ordered alignment object
    """
    first_index = 0 if alignment[0][0][0] == first else 1
    return [[item[first_index], item[1 ^ first_index]] for item in alignment]


def index2resi(protein):
    """returns a dictionary of [index] = resi"""
    res = {"res_i": dict()}
    cmd.iterate(protein, "res_i[index] = resv", space=res)
    return res["res_i"]


def atom_aln2resi_aln(raw_aln, target_ind=0):
    """
    Converts a raw_aln (parsed alignment object) to a dictionary mapping
    aligned residues between the target and the mobile protein.
    :param raw_aln:
    :param target_ind: 0 if target's atoms are the first tuple of every list in
    raw_aln, else 1
    :return: an ordered dictionary [ target residue ] = mobile residue
    """
    target_resi = index2resi(protein=raw_aln[0][target_ind][0])
    mobile_resi = index2resi(protein=raw_aln[0][1 ^ target_ind][0])
    resi_pairs = OrderedDict()
    for i in raw_aln:
        # X_r is the index's residue in protein X in the current position in raw_aln
        target_r = target_resi[i[target_ind][1]]
        mobile_r = mobile_resi[i[1 ^ target_ind][1]]
        if (target_r in resi_pairs.keys()) and (mobile_r != resi_pairs[target_r]):
            print(
                "conflict: target={:5<} mobile1={:5<} mobile2={:5<}".format(
                    target_r, resi_pairs[target_r], mobile_r
                )
            )
        else:
            resi_pairs[target_r] = mobile_r
    return resi_pairs


def get_alignment(mobile, target):
    """
    alignes mobile to target (both are loaded protein's names).
    :param mobile:
    :param target:
    :param cutoff:
    :return: the alignment object - list of lists of pairs of tuples:
    [(prot_a, index), (prot_b, index)]
    """
    cmd.select("ca_mobile", f"{mobile} and name CA")
    cmd.select("ca_target", f"{target} and name CA")
    aln_object = f"aln_{mobile}_{target}"
    cmd.cealign(mobile=mobile, target=target, object=aln_object)
    alignment = cmd.get_raw_alignment(aln_object)
    alignment = order_alignment(alignment, target)

    cmd.delete("ca_target")
    cmd.delete("ca_mobile")

    return alignment


def get_aligned_positions(template, mobile, positions: list):
    """"""
    template_name = load_protein(template)
    mobile_name = load_protein(mobile)

    alignment = get_alignment(mobile=mobile_name, target=template_name)
    res_aln = atom_aln2resi_aln(alignment, 0)

    return {int(pos): res_aln.get(int(pos)) for pos in positions}


class PymolAligner:
    def __init__(self, target: str, mobiles: list[str]) -> None:
        self.target = target
        self.mobiles = mobiles
        self.rmsds: dict[str:float] = {}

    def __call__(self) -> None:
        cmd.delete("all")
        cmd.load(self.target, "target")
        for ind in range(len(self.mobiles) // 500 + 1):
            from_ = ind * 500
            to_ = (ind + 1) * 500
            iter_names = self.mobiles[from_:to_]
            for i_name in iter_names:
                cmd.load(i_name, "mobile")
                self.rmsds[i_name] = cmd.align("mobile", "target")[0]
