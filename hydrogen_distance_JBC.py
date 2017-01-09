# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 18:22:34 2016

@author: hliu
"""

import os
import itertools
import mdtraj as md
import numpy as np


HOMEDIR = os.sep.join([os.environ['HOME'], 'Desktop', 'BBL_Salt_MakeUp'])
WORKDIR = os.sep.join([HOMEDIR, 'analysis', 'hydro_distance'])
os.chdir(WORKDIR)

traj = md.load('native.dcd', top=HOMEDIR+'/data/1w4h_native.pdb')

selected_residues = [14, # LEU
                     17, # HIS
                     19, # LEU
                     21, # ALA
                     24, # ILE
                     33, # LEU
                     38, # VAL
                     41, # HIS
                     42  # LEU
                    ]

region_seleRes = ['H1',
                  'H1',
                  'L',
                  '310',
                  '310',
                  'T',
                  'H2',
                  'H2',
                  'H2']

selected_resIDX = np.array(selected_residues) - 1



selected_atoms = [[221, 222, 223, 225, 226, 227], # HD
                  [265, 271],                     # HD
                  [298, 299, 300, 302, 303, 304], # HD
                  [325, 326, 327],                # HB
                  [365, 366, 367],                # HD
                  [484, 485, 486, 488, 489, 490], # HD
                  [567, 568, 569, 571, 572, 573], # HG
                  [623, 629],                     # HD
                  [642, 643, 644, 646, 647, 648]  # HD
                 ]

selected_atomIDX = [np.array(a)-1 for a in selected_atoms]

calculated = {x: y for x, y in zip(selected_resIDX, selected_atomIDX)}
regions = {x: y for x, y in zip(selected_resIDX, region_seleRes)}

results = {}
for res in itertools.combinations(selected_resIDX, 2):
    print regions[res[0]], regions[res[1]]
    x_atoms = calculated[res[0]]
    y_atoms = calculated[res[1]]
    jbc_IDx, jbc_IDy = res[0]-3, res[1]-3
    atom_pairs = [(i, j) for i in x_atoms for j in y_atoms]
    ave_dist = []
    for p in atom_pairs:
        dist = md.compute_distances(traj, [list(p)]).mean()
        ave_dist.append(dist)
    results[(jbc_IDx, jbc_IDy)] = np.mean(ave_dist)

