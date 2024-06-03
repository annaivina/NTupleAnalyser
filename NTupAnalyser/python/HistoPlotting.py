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


data = TFile('/afs/cern.ch/user/a/aivina/aivina/Hc/Analyser_Hc/run/HcSignal/hist-MGPy8Hc.root', 'READ')


#Histograms to be plotted
#h_pt = data.Get('chtruefid_pt')
#h_y  = data.Get('chtruefid_y')
#h_dr = data.Get('chtruefid_dr')

#gam1_pt = data.Get('cp1truefid_pt')
#gam2_pt = data.Get('cp2truefid_pt')
#gam1_eta = data.Get('cp1truefid_eta')
#gam2_eta = data.Get('cp2truefid_eta')


h_pt = data.Get('ptyy')
h_y  = data.Get('yyy')
h_m = data.Get('myy')

gam1_pt = data.Get('y1pt')
gam2_pt = data.Get('y2pt')
gam1_eta = data.Get('y1eta')
gam2_eta = data.Get('y2eta')


#pt
c1 = TCanvas('#gammagamma pT', '#gammagamma pT', 800, 600)
c1.cd()
c1.cd().SetLogy()
MakeColor(h_pt, 20, 4)
h_pt.GetXaxis().SetTitle("p_{#gamma #gamma}^{T} [GeV]")
h_pt.GetYaxis().SetTitle("Entries")
h_pt.Draw()
c1.SaveAs('gamgam_pt_reco.pdf')

#rapidity
c2 = TCanvas('#gammagamma y', '#gammagamma y', 800, 600)
c2.cd()
c2.cd().SetLogy()
MakeColor(h_y, 20, 4)
h_y.GetXaxis().SetTitle("y_{#gamma #gamma}")
h_y.GetYaxis().SetTitle("Entries")
h_y.Draw()
c2.SaveAs('gamgam_y_reco.pdf')

#dR
c3 = TCanvas('#gammagamma dR', '#gammagamma dR', 800, 600)
c3.cd()
MakeColor(h_m, 20, 4)
h_m.GetXaxis().SetTitle("m_{#gamma #gamma} [GeV]")
h_m.GetYaxis().SetTitle("Entries")
h_m.Draw()
c3.SaveAs('gamgam_m_reco.pdf')

#pts 2 gammas
c4 = TCanvas('gammas pTs', 'gammas pTs', 800, 600)
c4.cd()
c4.cd().SetLogy()
MakeColor(gam1_pt, 20, 4)
MakeColor(gam2_pt, 20, 2)
gam1_pt.GetXaxis().SetTitle("p_{T}^{gam,fid} [GeV]")
gam1_pt.GetYaxis().SetTitle("Entries")
gam1_pt.Draw()
gam2_pt.Draw('same ][')

leg = TLegend(.55,.65,.75,.85)
leg.AddEntry(gam1_pt, '  Leading photon', 'LP')
leg.AddEntry(gam2_pt, '  Subleading photon', 'LP')

leg.SetTextFont(42);
leg.SetTextSize(0.04);
leg.SetBorderSize(0);
leg.SetFillColor(10);
leg.Draw('same')
c4.SaveAs('gammas_pt_reco.pdf')

#etas 2 gammas
c5 = TCanvas('gammas etas', 'gammas etas', 800, 600)
c5.cd()
c5.cd().SetLogy()
MakeColor(gam1_eta, 20, 4)
MakeColor(gam2_eta, 20, 2)
gam1_eta.GetXaxis().SetTitle("#eta^{gam,fid}")
gam1_eta.GetYaxis().SetTitle("Entries")
gam1_eta.Draw()
gam2_eta.Draw('same ][')

c5.SaveAs('gammas_eta_reco.pdf')
