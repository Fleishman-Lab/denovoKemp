#!/usr/bin/env python3

import os
import logging
import argparse
import sys
import pymol
from pymol import cmd

__author__ = 'Dina Listov'


LGR = logging.getLogger(__file__)

"""usage struct_align.py <to be moved> <template>"""

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

def align(mobile, target):
    """All params are names of selections. The target will move on the PyMOL session
    """
    #score = cmd.align(mobile=mobile , target=target)
    score = cmd.cealign(target=target , mobile=mobile)
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

template_path = sys.argv[-1]
mobile_path = sys.argv[-2]

mobile_name= load_protein(mobile_path)
template_name=load_protein(template_path)
align(mobile_name, template_name)

index = mobile_path.rfind('.pdb')
save_path= os.getcwd()+'/'+ os.path.basename(mobile_path[:index])


save_structure(mobile_name, path=mobile_path)
