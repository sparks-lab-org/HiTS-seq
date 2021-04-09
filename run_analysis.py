#!/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
import argparse
import csv
import numpy as np
from collections import defaultdict
from sys import path
from hits.statistics.analysis import fit_background, analyze_samples
from hits.statistics.data import read_data
from hits.statistics.plots import position_plot
from hits.structure.secstruc import get_ss_blocks, read_fasta, read_ddomain
from hits.structure.secstruc_render import render_ss
from hits.structure.pymol import PymolData
import xlsxwriter as xl

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('gene', help='gene_name')
ap.add_argument('bdir', help='data directory')
ap.add_argument('fasta', type=str, default=None, help='fasta file')
ap.add_argument('dssp', type=str, default=None, help='dssp file')
ap.add_argument('ddomain', type=str, default=None, help='ddomain file')
ap.add_argument('--boundary', type=int, nargs='*')
ap.add_argument('--pymol', type=str, nargs='*')
ap.add_argument('--outdir', type=str, default='outputs/')

cmd = ap.parse_args()

gene = cmd.gene
bdir = cmd.bdir

########### READ DATA ########################
print("Reading Data: %s ..."%gene)
data = read_data(bdir, gene)
gene_length = len(data['A']['Antiparallel'][0])
nterm = 20
cterm = gene_length - 20
refseq = read_fasta(cmd.fasta)[gene]
##############################################

if cmd.pymol:
    if len(cmd.pymol) == 3:
        #pymol, pdb, chain
        pymol_dat = PymolData(cmd.pymol[0], cmd.pymol[1], refseq, chain=cmd.pymol[2])
    elif len(cmd.pymol) == 2:
        pymol_dat = PymolData(cmd.pymol[0], cmd.pymol[1], refseq)
    else:
        print("Invalid PYMOL")
        from sys import exit
        exit(1)

print("Starting analysis: %s ..."%gene)
fig, ax = plt.subplots(4, 1, figsize=(15,8), gridspec_kw={"height_ratios": [1,15,15,7.5], "hspace": 0.1})
#hack to avoid problems with sharex formatting
ax[0].set_xlim([0, gene_length+1])

for i in range(3):
    plt.setp(ax[i].get_xticklabels(), visible=False)
plt.setp(ax[3].get_yticklabels(), visible=False)

workbook = xl.Workbook(cmd.outdir+'/%s.xlsx'%gene)
bold = workbook.add_format({'bold': True})

#### Parameters #######
### dispersion: pseudo-replicates: All Antiparallel A/R (jointly with dilution)
sample1 = np.reshape(data['A']['Antiparallel'], -1)
sample2 = np.reshape(data['R']['Antiparallel'], -1) 
mask = (sample1>0) & (sample2>0)
theta0, beta0 = fit_background((sample1[mask], sample2[mask]))

#### TEST1 - SELF-COMPLEMENTATION (Parallel-A vs Parallel-T) ######
#Given inter-sample dispersion, re-estimate dilution factor between A/T library
### dilution: pseudo-replicates: ALL Antiparallel A/T (excluding termini)
sample1 = np.reshape(data['T']['Antiparallel'][:,nterm:cterm], -1)
sample2 = np.reshape(data['A']['Antiparallel'][:,nterm:cterm], -1) 
mask = (sample1>0) & (sample2>0)
theta, beta = fit_background((sample1[mask], sample2[mask]), theta=theta0)
#theta, beta = fit_background((sample1[mask], sample2[mask]))

print("Self:")
print("Beta = %.4f" % np.log2(beta))
print("Theta = %.4f" % np.log10(theta))

#Evaluate functional variants
sample1 = np.reshape(data['T']['Parallel'][0,:],-1)
sample2 = np.reshape(data['A']['Parallel'][0,:],-1)
log_fc, observed_mask, pv2, observed, unobserved = \
        analyze_samples(sample1, sample2, theta, beta)
#plot variants
position_plot(log_fc[observed_mask], pv2[observed_mask], \
        observed, unobserved, ax[1], 'Self')
ax[1].plot(range(1, gene_length+1), [np.log2(beta)]*gene_length, color='k', linestyle='--')

if cmd.pymol:
    pymol_dat.prepare_pymol(pv2, 'Self')

#Log results
worksheet = workbook.add_worksheet("Self")
header = ('Index', 'Count T - AP', 'Count T - P', 'LFC', 'log10(Pvalue)')
idx = [i+1 for i in range(len(log_fc))]
formatting = (None, None, None, "%.2f", "%.2f")

for col, entry in enumerate(header):
    worksheet.write(0, col, entry, bold)
for row, info in enumerate(sorted(zip(idx, sample1, sample2, log_fc, pv2), key = lambda x: x[-1])):
    for col, entry in enumerate(info):
        if formatting[col]:
            worksheet.write(row+1, col, float(formatting[col]%entry))
        else:
            worksheet.write(row+1, col, entry)

#### TEST2 - INDUCED-COMPLEMENTATION (Parallel-R vs Parallel-A) ######
#Global parameters already fit from A/R antiparallel
theta, beta = theta0, beta0

print("Assisted:")
print("Beta = %.4f" % np.log2(beta))
print("Theta = %.4f" % np.log10(theta))

#Evaluate functional variants
sample1 = np.reshape(data['A']['Parallel'][0,:],-1)
sample2 = np.reshape(data['R']['Parallel'][0,:],-1)
log_fc, observed_mask, pv2, observed, unobserved = \
        analyze_samples(sample1, sample2, theta, beta)

#Plot functional variants
position_plot(log_fc[observed_mask], pv2[observed_mask], \
        observed, unobserved, ax[2], 'Assisted')
ax[2].plot(range(1, gene_length+1), [np.log2(beta)]*gene_length, color='k', linestyle='--')
if cmd.pymol:
    pymol_dat.prepare_pymol(pv2, 'Assisted')

#Log results
worksheet = workbook.add_worksheet("Induced")
header = ('Index', 'Count A - P', 'Count R - P', 'LFC', 'log10(Pvalue)')
idx = [i+1 for i in range(len(log_fc))]
formatting = (None, None, None, "%.2f", "%.2f")

for col, entry in enumerate(header):
    worksheet.write(0, col, entry, bold)
for row, out in enumerate(sorted(zip(idx, sample1, sample2, log_fc, pv2), key = lambda x: x[-1])):
    for col, entry in enumerate(out):
        if formatting[col]:
            worksheet.write(row+1, col, float(formatting[col]%entry))
        else:
            worksheet.write(row+1, col, entry)
workbook.close()

### Add structure data to plots
ss_blocks, imap = get_ss_blocks(cmd.dssp, refseq)
render_ss(ss_blocks, ax[0])

idx, ddomain_dat = read_ddomain(cmd.ddomain, refseq, imap)
ax[3].plot(idx, ddomain_dat, color='k')
ax[3].set_ylabel('DDOMAIN')
ax[3].grid(axis='x', alpha=0.5)
ax[3].set_xlim([0, gene_length+1])
ax[3].xaxis.set_major_locator(ticker.MultipleLocator(20))
ax[3].xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax[3].get_yaxis().set_ticks([])
ymin, ymax = ax[3].get_ylim()
if cmd.boundary:
    for bnd0 in cmd.boundary:
        bnd = imap[bnd0]
        ax[3].plot([bnd, bnd], [ymin, ymax], linestyle='--', color='k')
fig.align_ylabels()

ax[1].grid(axis='x', alpha=0.5)
ax[2].grid(axis='x', alpha=0.5)
#fig.tight_layout()
plt.xlabel('Residue Position')
plt.savefig(cmd.outdir+'/'+gene+'_logfold_profile.png', bbox_inches='tight')

