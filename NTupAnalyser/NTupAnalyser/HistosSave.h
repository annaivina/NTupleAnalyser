#pragma once

/********************
 * HistogramStore:
 *   class to create, store, fill and retrieve TH1 and TH2 histograms
 *
 * Usage:
 *   HistogramStore HistoStore;
 *   HistoStore.createTH1F("Nphotons",40,-0.5,39.5,";#it{N}_{photon-clusters}");
 *   vector<TH1*> AllHistos = HistoStore.getListOfHistograms();
 *
 */

// STL include(s):
#include <map>
#include <vector>

// ROOT include(s):
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TString.h"


class HistosSave {
private:
  std::map<TString, TH1F *>     m_histoTH1F;
  std::map<TString, TH2F *>     m_histoTH2F;
  std::map<TString, TH3F *>     m_histoTH3F;
  std::map<TString, TProfile *> m_histoTProfile;

public:
  //create the histos
  void createTH1F(TString name, int Nbins, double xmin, double xmax, TString title = "");
  void createTH1F(TString name, const std::vector<double> &bins, TString title = "");
  void createTH2F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, TString title = "");
  void createTH2F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, TString title = "");
  void createTH3F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, int NBinsZ, double zmin, double zmax, TString title = "");
  void createTH3F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, const std::vector<double> &zbins, TString title = "");
  void createTProfile(TString name, int NbinsX, double xmin, double xmax, TString title = "");
  void createTProfile(TString name, const std::vector<double> &xbins, TString title = "");

  //filling the histos
  inline void fillTH1F(TString name, double x, double w = 1.0) {getTH1F(name)->Fill(x, w);}
  inline void fillTH2F(TString name, double x, double y, double w = 1.0) {getTH2F(name)->Fill(x, y, w);}
  inline void fillTH3F(TString name, double x, double y, double z, double w = 1.0) {getTH3F(name)->Fill(x, y, z, w);}
  inline void fillTProfile(TString name, double x, double y, double w = 1.0) {getTProfile(name)->Fill(x, y, w);}

  //checking the histos in the store
  inline bool hasTH1F(TString name) {return m_histoTH1F.count(name) > 0;}
  inline bool hasTH2F(TString name) {return m_histoTH2F.count(name) > 0;}
  inline bool hasTH3F(TString name) {return m_histoTH3F.count(name) > 0;}
  inline bool hasTProfile(TString name) {return m_histoTProfile.count(name) > 0;}

  //get the histos
  TH1F *getTH1F(TString name);
  TH2F *getTH2F(TString name);
  TH3F *getTH3F(TString name);

  //Make the profile
  TProfile *getTProfile(TString name);
  //getting the list of histograms
  std::vector<TH1 *> getListOfHistograms();
};
