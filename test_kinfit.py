#! /usr/bin/env python3
import numpy as np
from  KinFit import kinfit_3pr, fastmtt_cpp
import ROOT
from argparse import ArgumentParser

def fill_hist(hist, array):
    [hist.Fill(x) for x in array]

################
# Main routine #
################
parser = ArgumentParser()
parser.add_argument('-channel','--channel',dest='channel',default='tt',choices=['mt','tt'])

args = parser.parse_args()

# location of tuple
dirname='/eos/cms/store/group/phys_tau/lrussell/forAliaksei/CPSignalStudies/Run3_2022EE'
filename=dirname+'/'+args.channel+'/GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay/nominal/merged.root'

print('')
print('opening file %s'%(filename))    
df = ROOT.RDataFrame("ntuple",filename)
cuts = 'os>0.5&&idDeepTau2018v2p5VSe_2>=6&&idDeepTau2018v2p5VSmu_2>=4&&idDeepTau2018v2p5VSjet_2>=7&&pt_2>20.&&fabs(eta_2)<2.5&&decayModePNet_2==10&&hasRefitSV_2'

if args.channel=='mt':
    cuts += '&&iso_1<0.10&&pt_1>32&&fabs(eta_1)<2.1'
if args.channel=='tt':
    cuts += '&&idDeepTau2018v2p5VSe_1>=6&&idDeepTau2018v2p5VSmu_1>=4&&idDeepTau2018v2p5VSjet_1>=7&&pt_1>40.&&pt_2>40.&&fabs(eta_1)<2.5'
    
print('')
print('reading tuple as numpy columns')
cols = df.Filter(cuts).AsNumpy(["pt_1","pt_2","eta_1","eta_2","phi_1","phi_2","mass_1","mass_2","met_pt","met_phi","met_covXX","met_covXY","met_covYY","m_vis","sv_x_2","sv_y_2","sv_z_2","sv_cov00_2","sv_cov10_2","sv_cov11_2","sv_cov20_2","sv_cov21_2","sv_cov22_2","genPart_pt_1","genPart_eta_1","genPart_phi_1","genPart_pt_2","PVBS_x","PVBS_y","PVBS_z"])

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
results_MTT = fastmtt_cpp(cols["pt_2"],cols["eta_2"],cols["phi_2"],cols["mass_2"],cols['decay_type_1'],
                          cols["pt_1"],cols["eta_1"],cols["phi_1"],cols["mass_1"],cols['decay_type_2'],
                          cols["met_pt"],cols["met_phi"],
                          cols["met_covXX"],cols["met_covXY"],cols["met_covYY"],
                          mX,width)

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

####################################
#### Accessing results of fastMTT  #
####################################
x1 = results_MTT['x1'] # obtained w/o mass constraint
x2 = results_MTT['x2'] # obtained w/o mass constraint

x1_const = results_MTT['x1_BW'] # obtained with mass window cut
x2_const = results_MTT['x2_BW'] # obtained with mass window cut

pt1_fastMTT = cols['pt_1']/x1
pt2_fastMTT = cols['pt_2']/x2

pt1_fastMTT_const = cols['pt_1']/x1_const
pt2_fastMTT_const = cols['pt_2']/x2_const

mass = results_MTT['mass']

dpt1_kinfit = pt1/cols['genPart_pt_1']
dpt1_FastMTT = pt1_fastMTT/cols['genPart_pt_1'] 
dpt1_FastMTT_const = pt1_fastMTT_const/cols['genPart_pt_1']

dpt2_kinfit = pt2/cols['genPart_pt_2'] 
dpt2_FastMTT = pt2_fastMTT/cols['genPart_pt_2'] 
dpt2_FastMTT_const = pt1_fastMTT_const/cols['genPart_pt_2']

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
fill_hist(hist_mass_FastMTT,mass)

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

