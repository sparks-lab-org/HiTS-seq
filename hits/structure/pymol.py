#!/bin/env python

import numpy as np
from Bio.PDB import PDBParser
from Bio.Data.SCOPData import protein_letters_3to1 as p321
from .secstruc import indexmap
from ..vis_config import generate_colors, SIG_THRESH

p = PDBParser()

class PymolData(object):
    def __init__(self, pymol_fn, pdb_fn, refseq, chain='A'):
        self.pymol_fn = pymol_fn

        with open(pymol_fn, 'w') as f:
            f.write('load %s\n'%pdb_fn)
            f.write('bg_color white\n')
            f.write('set ray_shadows, 0\n')
            f.write('set depth_cue, 0\n')
            #f.write('set cartoon_fancy_helices, 0\n')
            #f.write('set cartoon_cylindrical_helices, 1\n')
            f.write('hide all\n')
            f.write('show cartoon, chain %s\n'%chain)
            f.write('zoom chain %s\n'%chain)
            #f.write('set_color mygrey, [%.3f, %.3f, %.3f]\n'%tuple(NONSIG_COLOR[0:3]))
            f.write('set_color mygrey, [0.7, 0.7, 0.7]\n')
            f.write('color mygrey\n')
            #set light_count, <integer>  

        structure = p.get_structure("", pdb_fn)
        residues = [r for r in structure[0][chain].get_residues() if r.id[0].strip()=='']
        pdbseq = "".join([p321.get(r.resname,"X") for r in residues])

        self.resids = [int(r.id[1]) for r in residues]
        self.imap = indexmap(refseq, pdbseq)
        self.color_dict = {}

    def prepare_pymol(self, stats_corrected, ylab):

        stat_colors = generate_colors(stats_corrected, ylab)
        #imap_inv = {v:k for k,v in self.imap.items()}

        with open(self.pymol_fn, 'a') as f:
            #REMOVE "EXTRAS"
            #for idx, rid in enumerate(self.resids):
            #    if idx not in imap_inv:
            #        f.write('hide cartoon, resi %s\n'%rid)

            for idx, (col, pv) in enumerate(zip(stat_colors, stats_corrected)):
                if pv>SIG_THRESH: continue
                if idx in self.imap:
                    key = tuple(col[0:3])
                    if key in self.color_dict:
                        colname = self.color_dict[key]
                    else:
                        colname = 'mycol'+str(len(self.color_dict.keys()))
                        self.color_dict[key] = colname
                        f.write('set_color %s, [%.3f, %.3f, %.3f]\n'%(colname,key[0], key[1], key[2]))
                    f.write('color %s, resi %d\n'%(self.color_dict[key], self.resids[self.imap[idx]]))