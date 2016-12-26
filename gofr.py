# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 20:39:03 2016

@author: hliu
"""

from plot_set import *
import glob
import os
import numpy as np
import scipy.io as sio
from matplotlib.ticker import FuncFormatter

concs = ['4M', '6M']
trajs = ['rg_large', 'rmsd_low']
labels = {'rg_large': 'Rg large', 'rmsd_low': 'RMSD low'}
ions = ['Li', 'Cl']
resnum = 45
cutoff = [1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 4.5, 5, 6]
# cutoff = [2.5]
cutoff_idx = lambda x: int(round(x/0.1))

def loadData(conc, traj, ion):
    path = os.environ['HOME']+os.sep+'Desktop/BBL_Salt_MakeUp/analysis/gofr/'+conc+os.sep+traj
    os.chdir(path)
    sarea = np.loadtxt('sarea_'+traj+'.dat')
    sarea = sarea.reshape((len(sarea)/resnum, resnum))
    unit_sarea = 1/sarea
    unit_sarea[np.isinf(unit_sarea)] = np.NaN
    salt_dist = sio.loadmat('ions_dist.mat')[ion+'_dist']
    return salt_dist, unit_sarea

fig_num = 0

y_formatter = FuncFormatter(lambda x, y: '%.1E' % x)

conc = '6M'
out1 = {}
out2 = {}

for ion in ions:
    for c in cutoff:
        fig_num += 1
        idx = cutoff_idx(c)
        f = plt.figure(fig_num)
        for t in trajs:
            salt_dist, unit_sarea = loadData(conc, t, ion)
            x = np.arange(1, 46)
            y = salt_dist[:, :, idx, 2]
            y = np.nanmean((y*unit_sarea), axis=0)
            if t == 'rg_large':
                out1[ion+'_C_'+str(c)] = y.reshape((len(y), 1))
            elif t == 'rmsd_low':
                out2[ion+'_C_'+str(c)] = y.reshape((len(y), 1))
            plt.plot(x, y, label=labels[t])
        plt.legend(frameon=False)
        plt.xlabel('Residue #')
        plt.ylabel(r'# of ions per $\AA^2$')
        plt.gca().yaxis.set_major_formatter(y_formatter)
        os.chdir(os.environ['HOME']+os.sep+'Desktop')
        plt.tight_layout()
        f.savefig('%s_%.2f.png' % (ion, c))
sio.savemat(conc+'_rg_large_'+'_ionDist_byCutoff.mat', out1)
sio.savemat(conc+'_rmsd_low_'+'_ionDist_byCutoff.mat', out2)