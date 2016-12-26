# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 18:22:34 2016

@author: hliu
"""

import itertools
import mdtraj as md

traj = md.load('protein.dcd', top='../1w4h_native.pdb')

selected_residues = [14,17,19,21,24,33,38,41,42]
selected_resIDX = np.array(selected_residues) - 1
for i in itertools.combinations(selected_resIDX, 2):
    print i

selected_atoms = [[221, 222, 223, 225, 226, 227],
                  [265, 271],
                  [298, 299, 300, 302, 303, 304],
                  [325, 326, 327],
                  [365, 366, 367],
                  [484, 485, 486, 488, 489, 490],
                  [567, 568, 569, 571, 572, 573],
                  [623, 629],
                  [642, 643, 644, 646, 647, 648]
                 ]

selected_atomIDX = [np.array(a)-1 for a in selected_atoms]

calculated = {x:y for x,y in zip(selected_resIDX, selected_atomIDX)}

results = {}
for res in itertools.combinations(selected_resIDX, 2):
    x_atoms = calculated[res[0]]
    y_atoms = calculated[res[1]]
    jbc_IDx, jbc_IDy = res[0]-3, res[1]-3
    atom_pairs = [(i, j) for i in x_atoms for j in y_atoms]
    ave_dist = []
    for p in atom_pairs:
        dist = md.compute_distances(traj, [list(p)]).mean()
        ave_dist.append(dist)
    results[(jbc_IDx, jbc_IDy)] = np.mean(ave_dist)