#ifndef __HcAnalyse_h
#define __HcAnalyse_h

#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include "EventLoop/Job.h"
#include "EventLoop/Worker.h"
#include "EventLoop/Algorithm.h"

#include "xAODRootAccess/tools/Message.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include "PathResolver/PathResolver.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODTracking/TrackParticleContainer.h"

#include <xAODJet/JetContainer.h>
#include <xAODEgamma/PhotonContainer.h>
#include <xAODTruth/TruthParticleContainer.h>


//Local includes
#include "NTupAnalyser/Common.h"
#include "NTupAnalyser/TruthHelper.h"
#include "NTupAnalyser/HistosSave.h"


using namespace std;

class HcAnalyse : public EL::Algorithm
{
    private:
            std::unique_ptr<TH1F>   m_cutflow_ext; //!
            HistosSave *m_hStore; //!

    		// Cut-flow - need to keep the same order!
    		//enum CutEnum {ALL = 0, PRESEL = 1, MATCH = 2, TIGHTID = 3, ISO = 4, RELPT = 5, MASS = 6, PHPT = 7, JET1 = 8, LEADJET2p5 = 9, LEADJET25 = 10, LEADJETNOBTAG = 11, LEADJETCTAG = 12, DRYJ = 13, HighPt =14};


    		// names of all cuts (do not includ "pass all")
    		//const std::vector<TString> s_cutDescs =
    		//{"All","Preselected","Trigger Match","tight ID","isolation","rel. #it{p}_{T} cuts","#it{m }_{#gamma#gamma} #in [105,160] GeV","#gamma p_{T}>25 GeV","AtLeast 1 jet", "lead j |#eta|<2.5", "lead j p_{T}>25 GeV","not B jets", "c tagged", "dR_y_j","HighPtcut"};

			enum CutEnum {ALL = 0, CJET = 1, PRESEL = 2, SEL = 3, PHPT = 4, JET1 = 5 };


    		// names of all cuts (do not includ "pass all")
    		const std::vector<TString> s_cutDescs =
    		{"All","Preselected","isPassed", "#gamma p_{T}>25 GeV","AtLeast 1 jet"};


  			/// value of cut that fail selection: PASSALL if all cuts passed
  			CutEnum CutFlowNew;
       std::map<std::string, TH1 *>  m_cutflows; //!

	public:

		TTree * tree; //!

		int run_n, lbn_n;//!
		int m_eventCounter;//!
		std::string _outputName;

		//static const int bins=41;//!
    	//float  count[bins]={-5.0, -4.75, -4.5, -4.25, -4.0, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75,
    	// 					-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5,
    	// 					4.75, 5.0};//!

		static const int bins=101;//!
        float  count[bins]={-5.0,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,
							-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
							0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
							3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};//!


    	static const int fbins=101;//!
    	float  fcount[fbins]={0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,
    						   0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,
    						   0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,
    						   0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,
    						   0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,
    						   0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00};//!

    	//The new scan (with few bins)
    	//static const int bins=21;//!
    	//float  count[bins]={-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};//!

		//static const int fbins=41;//!
    	//float  fcount[fbins]={0.00, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45,
    	//					  0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925,
    	//					  0.95, 0.975, 1.00};//!

        //static const int bins=41;//!
    	//float  count[bins]={-5.0, -4.75, -4.5, -4.25, -4.0, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75,
    	// 					-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5,
    	 //     				4.75, 5.0};//!

		//static const int bins=21;//!
    	//float  count[bins]={-5.0, -4.5,-4.0, -3.5,-3.0,-2.5,-2.0,-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,4.0, 4.5, 5.0};//!

    	//static const int fbins=21;//!
    	//float  fcount[fbins]={0.00, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,0.9, 0.95,1.00};//!


		//static const int bins=23;//!
    	//float  count[bins]={0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5};//!

    	//static const int fbins=23;//!
    	//float  fcount[fbins]={0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37};//!


		TH1D *masses_mv60[bins][fbins];//!
		TH1D *masses_mv60nonc[bins][fbins];//!
		TH1D *masses_mv70[bins][fbins];//!
		TH1D *masses_mv77[bins][fbins];//!
		TH1D *masses_mv85[bins][fbins];//!

		TH1D *masses_mv60_nodl1[bins][fbins];//!

		TH1D *masses_mv60_jvt[bins][fbins];//!
		TH1D *masses_mv70_jvt[bins][fbins];//!
		TH1D *masses_mv77_jvt[bins][fbins];//!
		TH1D *masses_mv85_jvt[bins][fbins];//!

		TH1D *masses_mv60w[bins][fbins];//!
		TH1D *masses_mv70w[bins][fbins];//!
		TH1D *masses_mv77w[bins][fbins];//!
		TH1D *masses_mv85w[bins][fbins];//!

        TH1D *dl1_60[fbins];//!
        TH1D *dl1_70[fbins];//!
        TH1D *dl1_77[fbins];//!
        TH1D *dl1_85[fbins];//!


        TH1D *dl1_60w[fbins];//!
        TH1D *dl1_70w[fbins];//!
        TH1D *dl1_77w[fbins];//!
        TH1D *dl1_85w[fbins];//!


		TH1D *N_recob[bins][fbins];//!
		TH1D *N_recoc[bins][fbins];//!
		TH1D *N_recol[bins][fbins];//!

		TH1D *N_tagb[bins][fbins];//!
		TH1D *N_tagc[bins][fbins];//!
		TH1D *N_tagl[bins][fbins];//!


		TH1D *eff_bmiss[bins][fbins];//!
		TH1D *eff_lmiss[bins][fbins];//!
		TH1D *eff_cmiss[bins][fbins];//!

		TH1D *eff_bmissw[bins][fbins];//!
		TH1D *eff_lmissw[bins][fbins];//!
		TH1D *eff_cmissw[bins][fbins];//!

		TH1D *eff_bmiss_nosf[bins][fbins];//!
		TH1D *eff_lmiss_nosf[bins][fbins];//!
		TH1D *eff_cmiss_nosf[bins][fbins];//!

		TH1D *eff_bmiss_init[bins][fbins];//!
		TH1D *eff_lmiss_init[bins][fbins];//!
		TH1D *eff_cmiss_init[bins][fbins];//!


		TH2D *j_ctrue;//!
		TH2D *j_btrue;//!
		TH2D *j_ltrue;//!

		TH2D *c_ctrue;//!
		TH2D *c_btrue;//!
		TH2D *c_ltrue;//!






		//Initialise the histos
		//photons
		//TH1F *infoEvent;//!
		TH1F *hist_out;//!


		float higgs_pt;//!
    	float higgs_y;//!
    	float higgs_m;//!
    	float higgs_e; //!

    	int jet_is_b;//!


		float topoetcone40_leading;//!
    	float topoetcone40_subleading;//!

    	float topoetcone20_leading;//!
    	float topoetcone20_subleading;//!

    	float ptcone20_leading;//!
  		float ptcone40_leading;//!

  		float ptcone20_subleading;//!
  		float ptcone40_subleading;//!

  		float  pt_leading;//!
    	float  pt_subleading;//!

    	float dl1_ctag_discriminant;//!

    	bool isPassedPreselection;//!
  		bool isPassedTriggerMatch ;//!
 		bool isPassedPID;//!
  		bool isPassedRelPtCuts;//!
  		bool isPassedMassCut;//!
  		bool isGoodJet;//!
  		bool isAtLeast1Jet;//!
  		bool isIsolated;//!
  		bool is_passed_primary;//!

  		bool is_btag_60;//!
    	bool is_btag_70;//!
    	bool is_btag_77;//!
    	bool is_btag_85;//!

    	int tagged_cnotb;//!
    	int jet_index;//!






  		bool is_Passed_FixCutTightCaloOnly;//!
  		bool is_Passed_FixCutTight;//!

  		int N_jets;//!
  		int N_1j0c;//!
        int N_1j1c;//!
    	int N_1exj0c;//!
        int N_1exj1c;//!

        float mass_gev;//!

        bool good_jet_20;//!
  		bool good_jet_15;//!
  		bool good_jet_25;//!
  		bool good_jet_30;//!
  		bool good_jet_35;//!


  		bool dr05_nomhj;//!
  		bool dr10_nomhj;//!
  		bool dr15_nomhj;//!
	    bool dr05_mhj150;//!
	    bool dr05_mhj155;//!
	    bool dr05_mhj160;//!
	    bool dr05_mhj165;//!
	    bool dr10_mhj150;//!
	    bool dr10_mhj155;//!
	    bool dr10_mhj160;//!
	    bool dr10_mhj165;//!
	    bool dr15_mhj150;//!
	    bool dr15_mhj155;//!
	    bool dr15_mhj160;//!
	    bool dr15_mhj165;//!



		//Initialized here, used by config
		//if this is from config - u do not use //!
		bool is_Data;//!
		bool isMC;//!
		bool _ispTcut;
		bool _do_Cutflow;
		bool _istruth;
		bool _isreco;
		string _reco_jet_collection;
		string _b_tagger;
		string _submitDir;
		string _mcname;
		bool _signal;
		float _jetpT;
    	float _jeteta;
    	float _cut;
    	float _frac;
    	bool _Dooptimization;
    	float _phpTcut;
    	float _dR_cut;



		//double invGeV = 0.001;//!
		//double GeV = 1000;//!





        float initial_weight;//!
        float total_weights;//!
        float pileup_weight;//!
    	float vertex_weight;//!
    	float scale;//!
    	float weight_analyse;//!
    	float weight_analyse20;//!
    	float weight_analyse30;//!
    	float weight_analyse35;//!
    	float weight_analyse50;//!
    	float tot_weight_nosf;//!

    	float weight_analyse_jvt;//!
    	float weight_analyse_jvt_ctagnom;//!
    	float weight_analyse_jvt_light;//!

    	float weight_ctag;//!
    	float weight_nonctag;//!


    	float int_weights;//!
        float jvt_weights;//!
        float mis_weights;//!
        float tot_weights;//!
        float mass_yy;//!
        float mass_yy_pv;//!
        float mass_yy_true;//!
        float xsectionbr;//!


     private:
      		TH1F *makeCutFlowHisto(int id, TString suffix = "");
      		TH1F *makeEventHisto(int id,TString suffix = "");
            std::map<int, TH1F *> m_cFlowHistos;
            std::map<int, TH1F *> m_EventInfo;
            int getSampleID = 1;

            TH1F *getCutFlowHisto(bool withDalitz = true)
  			{
    			int ID = getSampleID * (withDalitz ? -1 : 1);//if true set to -1 else 1

    			if (TH1F *h = m_cFlowHistos[ID]) { return h; }

    			m_cFlowHistos[ID] = makeCutFlowHisto(ID, withDalitz ? "Cutflow" : "_noDalitz");
    			return m_cFlowHistos[ID];
  			}

  			 TH1F *getEventHisto()
  			{
    			int ID = 1;
    			if (TH1F *h = m_EventInfo[ID]) { return h; }
    			m_EventInfo[ID] = makeEventHisto(ID,"EventInfo" );
    			return m_EventInfo[ID];
  			}

         	/// \brief fill the cut flow histograms
  	        void fillCutFlow(CutEnum cut);
  			// apply cut flow
  			// returns enum corresponding to failed cut - or PASSALL if all cuts passed
  			//CutEnum cutflow();


    public:

		HcAnalyse();

		/// \brief print a given cut flow histogram
  	    void printCutFlowHisto(TH1F *h, int Ndecimals = 0);
  	    void WeightsScales(const xAOD::EventInfo* HgameventInfo);
  	    double discriminant(double pc, double pb, double pu, double fraction);
  	    int identify_jets(const xAOD::JetContainer *jets, const xAOD::JetContainer *Trujets);
  	    void saveTxt(float a, float b, float c, float d, float e);
  	    double multiplyJvtWeights(const xAOD::JetContainer *jets);



		// these are the functions inherited from Algorithm
		virtual EL::StatusCode setupJob (EL::Job& job);
		virtual EL::StatusCode fileExecute ();
		virtual EL::StatusCode histInitialize ();
		virtual EL::StatusCode changeInput (bool firstFile);
		virtual EL::StatusCode initialize ();
		virtual EL::StatusCode execute ();
		virtual EL::StatusCode postExecute ();
		virtual EL::StatusCode finalize ();
		virtual EL::StatusCode histFinalize ();
		virtual EL::StatusCode createHistos();
        virtual EL::StatusCode doreco(const xAOD::PhotonContainer *photons,const xAOD::JetContainer *jets);
        virtual EL::StatusCode dotruth(const xAOD::TruthParticleContainer *truphotons,const xAOD::EventInfo* HgamTrutheventInfo);
        virtual EL::StatusCode checkJets(const xAOD::EventInfo* HgameventInfo,const xAOD::JetContainer *jets, const xAOD::PhotonContainer *photon, float weightjvt, float weightjvt30, float weightjvt50);

	protected:
	        inline virtual HistosSave  *HStore()  { return m_hStore; }

		// this is needed to distribute the algorithm to the workers
		ClassDef(HcAnalyse, 1);
};

#endif
