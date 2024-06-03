import glob, sys, os
import commands
from ROOT import *
from glob import *

import time
from math import sqrt

gROOT.LoadMacro("../macros/AtlasStyle.C")
gROOT.LoadMacro("../macros/AtlasUtils.C")
SetAtlasStyle()


def MakeColor(h, MarkerStyle, MarkerColor) :

    h.SetMarkerStyle(MarkerStyle)
    h.SetMarkerColor(MarkerColor)
    h.SetLineColor(MarkerColor)

#Obtaining the samples
#not full
ggf = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/ggF_ctag_full/hist-ggF.root', 'READ')
myy = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/myy_90_175_ctag/hist-myy_90_175.root', 'READ')
vbf = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/VBF_ctag/hist-VBF.root', 'READ')
wph = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/Wplus_ctag/hist-WplusH.root', 'READ')
wmh = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/Wminus_ctag/hist-WminusH.root', 'READ')
ggzh = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/ggZH_ctag/hist-ggZH.root', 'READ')
qqzh = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/qqZH_ctag/hist-qqZH.root', 'READ')
#signal
sig = TFile('/srv01/cgrp/users/annai/annai/Hc_Yukawa_Analysis/Analyser/run/HcSignal_ctag/hist-MGPy8Hc.root', 'READ')
#-------------------------------


#Histograms to be plotted
yy_ggf = ggf.Get('m_yy_j25')
yy_myy = myy.Get('m_yy_j25')
yy_vbf = vbf.Get('m_yy_j25')
yy_wph = wph.Get('m_yy_j25')
yy_wmh = wmh.Get('m_yy_j25')
yy_ggzh = ggzh.Get('m_yy_j25')
yy_qqzh = qqzh.Get('m_yy_j25')
yy_sig = sig.Get('m_yy_j25')
#-------------------------------


#Normalisation factors
lumi      = 36.1
norm_ggf  = (110.1 * lumi)/8.1556e+07
#norm_myy  = (51823 * lumi * 137.57)/1.04725e+08 (for fast)
norm_myy  = (51823 * lumi)/
norm_vbf  = (8.578 * lumi)/3292930
norm_wph  = (1.902 * lumi)/200657
norm_wmh  = (1.206 * lumi)/125980
norm_qqzh = (1.725 * lumi)/339286
norm_ggzh = (0.2782* lumi)/2562.09
norm_sig  = (1.0192* lumi)/100000
#-------------------------------

#Scaling
yy_ggf.Scale(norm_ggf)
yy_myy.Scale(norm_myy)
yy_vbf.Scale(norm_vbf)
yy_wph.Scale(norm_wph)
yy_wmh.Scale(norm_wmh)
yy_ggzh.Scale(norm_ggzh)
yy_qqzh.Scale(norm_qqzh)
yy_sig.Scale(norm_sig)

print("The integral for ggF ")
print yy_ggf.Integral()

print("The integral for signal ")
print yy_sig.Integral()

print("The integral for myy ")
print yy_myy.Integral()

print("The integral for VBF")
print yy_vbf.Integral()

print("The integral for WpH ")
print yy_wph.Integral()

print("The integral for WmH ")
print yy_wmh.Integral()

print("The integral for qqZH ")
print yy_qqzh.Integral()

print("The integral for ggZH ")
print yy_ggzh.Integral()

#-------------------------------

#Plot all of them on one canvas
c1 = TCanvas('#gammagamma m', '#gammagamma m', 800, 600)
c1.cd()
c1.cd().SetLogy()
MakeColor(yy_ggf,  20, kRed)
MakeColor(yy_myy,  20, kOrange)
MakeColor(yy_vbf,  20, kYellow+2)
MakeColor(yy_wph,  20, kGreen+3)
MakeColor(yy_wmh,  20, kCyan+1)
MakeColor(yy_ggzh, 20, kAzure+4)
MakeColor(yy_qqzh, 20, kViolet-6)
MakeColor(yy_sig,  20, kBlack)
yy_ggf.GetXaxis().SetTitle("m_{#gamma #gamma} [GeV]")
yy_ggf.GetYaxis().SetTitle("Entries")
#Drawing
yy_ggf.Draw()
#yy_ggf.GetYaxis().SetRangeUser(e-04,e06)
yy_myy.Draw('same')
yy_vbf.Draw('same')
yy_wph.Draw('same')
yy_wmh.Draw('same')
yy_qqzh.Draw('same')
yy_sig.Draw('same')

#Legend
leg = TLegend(.55,.65,.75,.85)
leg.AddEntry(yy_myy, 'Sherpa #gamma#gamma', 'LP')
leg.AddEntry(yy_ggf, 'ggF', 'LP')
leg.AddEntry(yy_vbf, 'VBF', 'LP')
leg.AddEntry(yy_wph, 'W^{+}H', 'LP')
leg.AddEntry(yy_qqzh, 'qq#rightarrow ZH', 'LP')
leg.AddEntry(yy_wmh, 'W^{-}H', 'LP')
leg.AddEntry(yy_ggzh, 'gg#rightarrow ZH', 'LP')
leg.AddEntry(yy_sig, 'Signal', 'LP')






leg.SetTextFont(42);
leg.SetTextSize(0.04);
leg.SetBorderSize(0);
leg.SetFillColor(10);
leg.Draw('same')
c1.SaveAs('mass_jet_25.pdf')



#-------------------------------
