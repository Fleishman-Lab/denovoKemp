#!/bin/env python3

import argparse
import pymol_alignment
import os
import pymol
from pymol import cmd


def parse_args():
    """pymol script that returns aligned positions on the mutant structure given a template."""
    desc = ''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('template', help='path to a template\'s pdb file')
    parser.add_argument('mutant', help='path to a mutant\'s pdb file')
    parser.add_argument('positions', nargs='+', action='store',
                        help='space separated list of residue numbers')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    pos_dict = pymol_alignment.get_aligned_positions(template=args.template,
                                                       mobile=args.mutant,
                                                       positions=args.positions)
    output=' '.join('{}'.format(val) for val in pos_dict.values())
    print(output)

