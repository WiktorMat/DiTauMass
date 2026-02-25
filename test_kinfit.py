#! /usr/bin/env python3
import numpy as np
from  KinFit import kinfit_3pr#, fastmtt_cpp
import ROOT
from argparse import ArgumentParser

def fill_hist(hist, array):
    [hist.Fill(x) for x in array]

################
# Main routine #
################
parser = ArgumentParser()
parser.add_argument('input',help='input file (ROOT or CSV)')
parser.add_argument('-channel','--channel',dest='channel',default='tt',choices=['mt','tt'])

args = parser.parse_args()

# location of tuple
# dirname='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/CPSignalStudies/Run3_2022EE'
# filename=dirname+'/'+args.channel+'/GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay/nominal/merged.root'
filename = args.input
# filename = '/eos/home-w/wmatyszk/HiggsDNA/CleanDNA/higgs-dna-waw/combine_test_run/data/nominal.root'

print('')
print('opening file %s'%(filename))    
df = ROOT.RDataFrame("ntuple",filename)
cuts = 'os>0.5&&idDeepTau2018v2p5VSe_2>=6&&idDeepTau2018v2p5VSmu_2>=4&&idDeepTau2018v2p5VSjet_2>=7&&pt_2>20.&&fabs(eta_2)<2.5&&decayModePNet_2==10&&hasRefitSV_2'

if args.channel=='mt':
    cuts += '&&iso_1<0.10&&pt_1>32&&fabs(eta_1)<2.1'
if args.channel=='tt':
    cuts += '&&idDeepTau2018v2p5VSe_1>=6&&idDeepTau2018v2p5VSmu_1>=4&&idDeepTau2018v2p5VSjet_1>=7&&pt_1>40.&&pt_2>40.&&fabs(eta_1)<2.5'
    
print('')
# get list of all branches in the tree
all_cols = [str(c) for c in df.GetColumnNames()]

print(f'Found {len(all_cols)} columns')

print('')
print('reading tuple as numpy columns')

cols = df.Filter(cuts).AsNumpy(all_cols)

print('')
print('Length of column : %1i\n'%(len(cols["pt_1"])))
print('Running KinFit (be patient, it takes awhile)')

# steering parameters
phiScan = False # don't perform phi scan
mX = 125.10 # Higgs mass
width = 2.5 # Higgs mass window for constraint 

svx = cols['sv_x_2']-cols['PVBS_x']
svy = cols['sv_y_2']-cols['PVBS_y']
svz = cols['sv_z_2']-cols['PVBS_z']

cols['decay_type_1'] = 1*np.ones(len(cols["pt_1"]),dtype=np.int32)
cols['decay_type_2'] = 2*np.ones(len(cols["pt_1"]),dtype=np.int32)
if args.channel=='tt':
    cols['decay_type_1'] = 2*np.ones(len(cols["pt_1"]),dtype=np.int32)

#######################
# Calling kinfit_3pr  #
#######################
# The second tau is assumed to decay to a1 (see definition of cuts above!)
# should be passed first to the routine
results = kinfit_3pr(cols["pt_2"],cols["eta_2"],cols["phi_2"],cols["mass_2"],
                     cols["pt_1"],cols["eta_1"],cols["phi_1"],cols["mass_1"],
                     cols["met_pt"],cols["met_phi"],
                     cols["met_covXX"],cols["met_covXY"],cols["met_covYY"],
                     svx,svy,svz,
                     cols["sv_cov00_2"],cols["sv_cov10_2"],cols["sv_cov20_2"],
                     cols["sv_cov11_2"],cols["sv_cov21_2"],cols["sv_cov22_2"],
                     phiScan,mX)    
###################
# calling fastMTT #
###################
# results_MTT = fastmtt_cpp(cols["pt_2"],cols["eta_2"],cols["phi_2"],cols["mass_2"],cols['decay_type_1'],
#                           cols["pt_1"],cols["eta_1"],cols["phi_1"],cols["mass_1"],cols['decay_type_2'],
#                           cols["met_pt"],cols["met_phi"],
#                           cols["met_covXX"],cols["met_covXY"],cols["met_covYY"],
#                           mX,width)

###############################################################
# accessing results of kinfit_3pr (library: keyword->column)
###############################################################
# this is tau not decaying to a1 (second tau passed to kinfit)
px1 = results['px_2']
py1 = results['py_2']
pz1 = results['pz_2']

# this is tau decaying to a1 (first tau passed to kinfit)
px2 = results['px_1']
py2 = results['py_1']
pz2 = results['pz_1']

# chi-squared of the fit
chi2 = results['chi2']

pt1 = np.sqrt(px1**2+py1**2)
pt2 = np.sqrt(px2**2+py2**2)

# --- store KinFit outputs ---
cols["KinFit_pt_1"] = pt1
cols["KinFit_pt_2"] = pt2

####################################
#### Accessing results of fastMTT  #
####################################
# x1 = results_MTT['x1'] # obtained w/o mass constraint
# x2 = results_MTT['x2'] # obtained w/o mass constraint

# x1_const = results_MTT['x1_BW'] # obtained with mass window cut
# x2_const = results_MTT['x2_BW'] # obtained with mass window cut

# pt1_fastMTT = cols['pt_1']/x1
# pt2_fastMTT = cols['pt_2']/x2

# pt1_fastMTT_const = cols['pt_1']/x1_const
# pt2_fastMTT_const = cols['pt_2']/x2_const

# mass = results_MTT['mass']

dpt1_kinfit = pt1/cols['genPart_pt_1']
dpt1_FastMTT = cols['FastMTT_pt_1']/cols['genPart_pt_1'] 
dpt1_FastMTT_const = cols['FastMTT_pt_1_constraint']/cols['genPart_pt_1']

dpt2_kinfit = pt2/cols['genPart_pt_2'] 
dpt2_FastMTT = cols['FastMTT_pt_2']/cols['genPart_pt_2'] 
dpt2_FastMTT_const = cols['FastMTT_pt_2_constraint']/cols['genPart_pt_2']

# saving to RooT file
outputFile="kinfit_3pr_%s.root"%(args.channel)
f = ROOT.TFile(outputFile,"recreate")
f.cd('')

hist_dpt1_kinfit = ROOT.TH1D("dpt1_kinfit","",60,0.,3.)
hist_dpt2_kinfit = ROOT.TH1D("dpt2_kinfit","",60,0.,3.)
hist_dpt1_FastMTT = ROOT.TH1D("dpt1_fastMTT","",60,0.,3.)
hist_dpt2_FastMTT = ROOT.TH1D("dpt2_fastMTT","",60,0.,3.)
hist_dpt1_FastMTT_const = ROOT.TH1D("dpt1_fastMTT_const","",60,0.,3.)
hist_dpt2_FastMTT_const = ROOT.TH1D("dpt2_fastMTT_const","",60,0.,3.)
hist_mass_FastMTT = ROOT.TH1D('mass_FastMTT','',60,0.,300.)
hist_mvis = ROOT.TH1D('mvis','',60,0.,300.)

#####################
# Filling histograms
#####################

fill_hist(hist_dpt1_kinfit,dpt1_kinfit)
fill_hist(hist_dpt2_kinfit,dpt2_kinfit)
fill_hist(hist_dpt1_FastMTT,dpt1_FastMTT)
fill_hist(hist_dpt2_FastMTT,dpt2_FastMTT)
fill_hist(hist_dpt1_FastMTT_const,dpt1_FastMTT_const)
fill_hist(hist_dpt2_FastMTT_const,dpt2_FastMTT_const)
fill_hist(hist_mvis,cols['m_vis'])
# fill_hist(hist_mass_FastMTT,mass)

#####################
# saving histograms
#####################
f.cd('')
hist_mvis.Write("mvis")
hist_mass_FastMTT.Write("mtt")

hist_dpt1_kinfit.Write("dpt1_kinfit")
hist_dpt1_FastMTT.Write("dpt1_FastMTT")
hist_dpt1_FastMTT_const.Write("dpt1_FastMTT_const")

hist_dpt2_kinfit.Write("dpt2_kinfit")
hist_dpt2_FastMTT.Write("dpt2_FastMTT")
hist_dpt2_FastMTT_const.Write("dpt2_FastMTT_const")

f.Close()
print('')
print('Histograms are saved in file %s'%(outputFile))
print('')

################################
# Save new ROOT file with tree #
################################

print("\nSaving ROOT file with new columns...")

from pathlib import Path

# --- output path ---
outdir = Path.cwd()
outfile = outdir / "merged.root"

# 🔒 nie nadpisuj
if outfile.exists():
    i = 1
    while (outdir / f"merged_{i}.root").exists():
        i += 1
    outfile = outdir / f"merged_{i}.root"

########################################
# Create ROOT file and TTree manually #
########################################

import array

fout = ROOT.TFile(str(outfile), "RECREATE")
tree = ROOT.TTree("ntuple", "ntuple")

# --- prepare branch buffers ---
buffers = {}
branches = {}

skipped_cols = []

for name, arr in cols.items():
    dtype = arr.dtype

    if np.issubdtype(dtype, np.floating):
        buffers[name] = array.array('f', [0.])
        tree.Branch(name, buffers[name], f"{name}/F")
    elif np.issubdtype(dtype, np.integer):
        buffers[name] = array.array('i', [0])
        tree.Branch(name, buffers[name], f"{name}/I")
    else:
        # try to coerce object columns
        try:
            cols[name] = np.asarray(cols[name], dtype=np.int32)
            buffers[name] = array.array('i', [0])
            tree.Branch(name, buffers[name], f"{name}/I")
            print(f"[WARN] coerced {name} to int32")
        except Exception:
            try:
                cols[name] = np.asarray(cols[name], dtype=np.float32)
                buffers[name] = array.array('f', [0.])
                tree.Branch(name, buffers[name], f"{name}/F")
                print(f"[WARN] coerced {name} to float32")
            except Exception:
                skipped_cols.append(name)
                print(f"[SKIP] dropping non-numeric column: {name}")
                continue

n = len(next(iter(cols.values())))

print(f"Writing {n} entries...")

# --- event loop ---
active_cols = list(buffers.keys())

for i in range(n):
    for name in active_cols:
        buffers[name][0] = cols[name][i]
    tree.Fill()

print(f"Skipped {len(skipped_cols)} non-numeric columns: {skipped_cols}")

tree.Write()
fout.Close()

print(f"✅ ROOT saved to: {outfile}")
print("")