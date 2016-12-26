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

os.chdir('/home/hliu/Desktop/BBL_Salt_MakeUp/analysis/chemical_shift')

seq = np.loadtxt('sequence_3L.dat', dtype=str)
cs_native = sio.loadmat('native_chemicalShift.mat')
cs_random = np.loadtxt('cs_randomCoil.dat')
selected_atomTypes = ['H', 'N']

class ChemicalShift:

    def __init__(self, trjname, seq):
        self.name = trjname
        self.seq = pd.Series(seq,
                            index=np.arange(1, len(seq)+1))

    def loadShift(self, shiftData):
        atom_types = ['H', 'HA', 'N', 'C', 'CA', 'CB']
        shift_values = {t: shiftData[t].mean(axis=0) for t in atom_types}
        shift_values.update({t+'_E': shiftData[t].std(axis=0) for t in atom_types})
        shift_values = pd.DataFrame(shift_values)
        shift_values.index = np.arange(1, len(seq)+1)
        shift_values.index.name = 'ResID'
        self.shift = shift_values

    def setCoor(self, selected_atomTypes):
        self.coor = self.shift[selected_atomTypes]

    def setPlot(self, getPlot):
        self.plot = getPlot(self.coor)

plotted_index = np.delete(np.arange(1, 46), 7)
getPlot = lambda x: x.ix[plotted_idx, :]

def setData(trjName, fPth):
    cs = ChemicalShift(trjName, seq)
    cs.loadShift(sio.loadmat(fPth))
    cs.setCoor(selected_atomTypes)
    cs.setPlot(getPlot)
    return cs
s
cs_6 = setData('6M low', '6M/chemicalShift.mat')
cs_4 = setData('4M low', '4M/chemicalShift.mat')
cs_native = setData('Native', 'native_chemicalShift.mat')


comp1 = cs_4
comp2 = cs_native
plt.scatter(comp1.plot['H'], comp1.plot['N'], c='r')
plt.scatter(comp2.plot['H'], comp2.plot['N'], c='b')

for idx in plotted_idx:
    diff = np.linalg.norm(comp2.coor.loc[idx, :].values-comp1.coor.loc[idx, :].values)
    print diff, idx, comp1.seq[idx]

width = 0.3
comp1 = cs_4_coor
comp2 = cs_6_coor

plt.bar(comp1.index-width,
        comp1['HA'].values-cs_random,
        width=width,
        label=comp1.name)
plt.errorbar(comp1.index-width/2,
             comp1['HA'].values-cs_random,
             yerr=comp1['HA_E'].values,
             ecolor='gray',
             elinewidth=1.5,
             linewidth=0)
plt.bar(comp2.index,
        comp2['HA'].values-cs_random,
        width=width,
        color='r',
        label=comp2.name)
ax = plt.gca()
set_tickinterval(ax, 'x', 5)
plt.legend(loc=2, frameon=False)
plt.grid()
plt.xlim([-1,46])
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, i: str(int(x))))