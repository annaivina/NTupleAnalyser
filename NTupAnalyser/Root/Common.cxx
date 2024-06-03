// STL include(s):
#include <iostream>

// ROOT include(s):
#include "TFile.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include <TTree.h>
#include <vector>
#include <string>

// Local include(s):
#include "NTupAnalyser/Common.h"


// See header file for documentation
using namespace std;
namespace HC {


  	void fatal(TString msg){
    	printf("\nFATAL\n  %s\n\n", msg.Data());
    	abort();
  	}


 	void printCutFlowHisto(TH1F *h, int Ndecimals){
  		TString format("  %-24s%10.");
  		format += Ndecimals;
  		format += "f%11.2f%%%11.2f%%\n";
  		int all_bin = h->FindBin(0);
  		printf("  %-24s%10s%12s%12s\n", "Event selection", "Nevents", "Cut rej.", "Tot. eff.");

  		for (int bin = 1; bin <= h->GetNbinsX(); ++bin){
    		double ntot = h->GetBinContent(all_bin), n = h->GetBinContent(bin), nprev = h->GetBinContent(bin - 1);
    		TString cutName(h->GetXaxis()->GetBinLabel(bin));
    		cutName.ReplaceAll("#it{m}_{#gamma#gamma}", "m_yy");

    		if (bin == 1 || nprev == 0 || n == nprev)
    		{ printf(format.Data(), cutName.Data(), n, -1e-10, n / ntot * 100); }
    		else // if the cut does something, print more information
    		{ printf(format.Data(), cutName.Data(), n, (n - nprev) / nprev * 100, n / ntot * 100); }
  		}
	}

  TString fourVecAsText(const TLorentzVector &p4)
  {
    return TString::Format("(pT,y,phi,m) = (%6.1f GeV,%6.3f,%6.3f,%5.1f GeV )",
                           p4.Pt() * invGeV, p4.Rapidity(), p4.Phi(), p4.M() * invGeV);
  }

  TString fourVecAsText(const xAOD::IParticle *p)
  {
    return fourVecAsText(p->p4());
  }



  StrV vectorize(TString str, TString sep){
    StrV result;
    TObjArray *strings = str.Tokenize(sep.Data());

    if (strings->GetEntries() == 0) { delete strings; return result; }

    TIter istr(strings);

    while (TObjString *os = (TObjString *)istr()){
      // the number sign and everything after is treated as a comment
      if (os->GetString()[0] == '#') { break; }
      result.push_back(os->GetString());
    }
    delete strings;
    return result;
  }


  	// checks if a given file or directory exist
  	bool fileExist(TString fn){
    	return !(gSystem->AccessPathName(fn.Data()));
  	}


	TH1 *getHistogramFromFile(TString fname, TString hname){
    	fname = PathResolverFindCalibFile(fname.Data());
    	TFile *file = TFile::Open(fname.Data(), "READ");

    	if (file == nullptr){
      		std::cout << "HgammaUtils::getHistogramFromFile() : Couldn't open file "
            << fname.Data() << ", returning nullptr." << std::endl;
      		return nullptr;
    	}

    	TH1 *temp = dynamic_cast<TH1 *>(file->Get(hname.Data()));

    	if (temp == nullptr){
			std::cout << "HgammaUtils::getHistogramFromFile() : Couldn't find histogram "
    		<< hname.Data() << " in file "
			<< fname.Data() << ", returning nullptr." << std::endl;
      	return nullptr;
    	}

    	bool status = TH1::AddDirectoryStatus();
    	TH1::AddDirectory(false);
    	hname = "cloned_" + hname;
    	TH1 *hist = dynamic_cast<TH1 *>(temp->Clone(hname.Data()));
    	SafeDelete(file);
    	TH1::AddDirectory(status);

    	return hist;
	}
}
