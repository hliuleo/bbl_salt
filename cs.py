# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:31:02 2016

@author: hliu
"""

import scipy.io as sio
import pandas as pd
from matplotlib.ticker import FuncFormatter
from plot_set import *
import os

HOME = os.environ['HOME']
WORKDIR = HOME+os.sep+'Desktop'
CS_PTH = 'BBL_Salt_MakeUp/analysis/chemical_shift'
METHOD = 'ShiftX'
os.chdir(os.sep.join([WORKDIR, CS_PTH, METHOD]))

seq_3L = np.loadtxt('sequence_3L.dat', dtype=str)
seq_L = np.loadtxt('sequence.dat', dtype=str)
seqs = {'1L': seq_L, '3L': seq_3L}
num_res = len(seq_L)
cs_native = sio.loadmat('native_chemicalShift.mat')
cs_random = np.loadtxt('cs_randomCoil.dat')


class ChemicalShift:

    def __init__(self, trjname, seqs):
        self.name = trjname
        self.seq = pd.DataFrame(seqs, 
                                index=np.arange(1, num_res+1))

    def loadShift(self, shiftData):

        atom_types = ['H', 'HA', 'N', 'C', 'CA', 'CB']
        shift_values = {t: np.nanmean(shiftData[t], axis=0) for t in atom_types}
        shift_values.update({t+'_E': np.nanstd(shiftData[t], axis=0) for t in atom_types})
        shift_values = pd.DataFrame(shift_values)
        shift_values.index = np.arange(1, num_res+1)
        shift_values.index.name = 'ResID'
        self.shift = shift_values

    def setCoor(self, selected_atomTypes):
        self.coor = self.shift[selected_atomTypes]

    def setPlot(self, plotted_idx=None):
        if plotted_idx is None:
            self.plot = self.coor
        else:
            self.plot = self.coor.ix[plotted_idx, :]

plotted_idx = np.delete(np.arange(1, 46), 7)

def setData(trjName, fPth):
    cs = ChemicalShift(trjName, seqs)
    cs.loadShift(sio.loadmat(fPth))
    cs.setCoor(selected_atomTypes)
    cs.setPlot(plotted_idx)
    return cs

selected_atomTypes = ['H', 'N']
cs_6 = setData('6M low', '6M/chemicalShift.mat')
cs_4 = setData('4M low', '4M/chemicalShift.mat')
cs_native = setData('Native', 'native_chemicalShift.mat')


comp1 = cs_6
comp2 = cs_native
plt.scatter(comp1.plot['H'], comp1.plot['N'],
            c='b', alpha=0.5, label=comp1.name, s=100)
plt.scatter(comp2.plot['H'], comp2.plot['N'],
            c='r', alpha=0.5, label=comp2.name, s=100)

ax = plt.gca()
ax.invert_yaxis()
ax.invert_xaxis()

plt.legend(frameon=False, loc=2)

plotted_residue = [9, 17, 18, 20, 21, 24, 26, 30, 31, 34, 39, 40] # numbers are coming from and compatible with JBC

for res in plotted_residue:
    v1 = comp1.coor.ix[res+4].values
    v2 = comp2.coor.ix[res+4].values
    d = np.linalg.norm(v1-v2)
    plt.plot([v1[0], v2[0]],
             [v1[1], v2[1]],
             c='k', alpha=0.5)
    plt.text(v2[0], v2[1], comp1.seq['1L'][res+4]+str(res), fontsize=20)

'''
Draw lines for residues having large deviations

for idx in comp1.coor.index:
    v1 = comp1.coor.ix[idx].values
    v2 = comp2.coor.ix[idx].values
    d = np.linalg.norm(v1-v2)
    if d > 3.5:
        plt.plot([v1[0], v2[0]],
                 [v1[1], v2[1]],
                 c='k', alpha=0.5)
        plt.text(v2[0], v2[1], comp1.seq['1L'][idx]+str(idx-4), fontsize=20)
'''

plt.xlabel(r'$^{1}$H Chemical Shift')
plt.ylabel(r'$^{15}$N Chemical Shift')

for idx in plotted_idx:
    diff = np.linalg.norm(comp2.coor.loc[idx, :].values-comp1.coor.loc[idx, :].values)
    print diff, idx, comp1.seq['3L'][idx]

width = 0.3

selected_atomTypes = ['HA', 'HA_E']
cs_6.setCoor(selected_atomTypes)
cs_4.setCoor(selected_atomTypes)
cs_native.setCoor(selected_atomTypes)
comp1 = cs_4
comp2 = cs_6

plt.bar(comp1.coor.index-width,
        comp1.coor['HA'].values-cs_random,
        width=width,
        label=comp1.name)
plt.errorbar(comp1.coor.index-width/2,
             comp1.coor['HA'].values-cs_random,
             yerr=comp1.coor['HA_E'].values,
             ecolor='gray',
             elinewidth=1.5,
             linewidth=0)
plt.bar(comp2.coor.index,
        comp2.coor['HA'].values-cs_random,
        width=width,
        color='r',
        label=comp2.name)
ax = plt.gca()
set_tickinterval(ax, 'x', 5)
plt.legend(loc=2, frameon=False)
plt.grid()
plt.xlim([-1,46])
plt.xlabel('Residue ID')
plt.ylabel('Secondary Shift')
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, i: str(int(x))))
