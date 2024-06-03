//#include <AsgTools/MessageCheck.h>
//Local includes
#include "NTupAnalyser/Common.h"
#include "NTupAnalyser/TruthHelper.h"
#include "NTupAnalyser/HistosSave.h"
#include "NTupAnalyser/HcAnalyse.h"

//#include <AsgTools/MessageCheck.h>
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODEgamma/EgammaContainer.h>
#include <xAODEgamma/ElectronContainer.h>
#include <xAODEgamma/PhotonContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODEgamma/EgammaEnums.h>
#include <xAODPrimitives/IsolationType.h>
#include <xAODTracking/TrackingPrimitives.h>
#include <xAODTracking/TrackParticlexAODHelpers.h>
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"

// Infrastructure includes:
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/ShallowAuxContainer.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"

//Root includes
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <vector>
#include <string>



// this is needed to distribute the algorithm to the workers
ClassImp(HcAnalyse)


#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
  if( ! EXP.isSuccess() ) {				\
    Error( CONTEXT,					\
    XAOD_MESSAGE( "Failed to execute: %s" ),	\
    #EXP );					\
    return EL::StatusCode::FAILURE;			\
    }							\
    } while( false )



//Constructor
HcAnalyse :: HcAnalyse (void)
: m_hStore(nullptr)
{ }

EL::StatusCode HcAnalyse :: setupJob (EL::Job& job)
{
  job.useXAOD ();
  EL_RETURN_CHECK( "setupJOB()", xAOD::Init());
  return EL::StatusCode::SUCCESS;
}

//Initialising the histos for the Plots
EL::StatusCode HcAnalyse :: histInitialize ()
{
	//Initializing the histogram store
	m_hStore = new HistosSave();

	return EL::StatusCode::SUCCESS;
}


EL::StatusCode HcAnalyse :: fileExecute ()
{

   // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  TFile *output = dynamic_cast<TFile *>(wk()->getOutputFile(_outputName));

  // Copy cut flow and other histograms from input file
  TFile *input = dynamic_cast<TFile *>(wk()->inputFile());
  TIter next(input->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)next())) {
    std::string className = key->GetClassName();

    if (className.find("TH1") != std::string::npos) {
      const char *name = key->GetName();
      TH1 *cutflow = dynamic_cast<TH1 *>(wk()->inputFile()->Get(name));

      if (m_cutflows.find(name) == m_cutflows.end()) {
        m_cutflows[name] = dynamic_cast<TH1 *>(cutflow->Clone());
        m_cutflows[name]->SetDirectory(output);
      } else {
        m_cutflows[name]->Add(cutflow);
      }
    }
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HcAnalyse :: changeInput (bool /*firstFile*/)
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HcAnalyse :: initialize ()
{
	m_eventCounter = 0;
	xAOD::TEvent* event = wk()->xaodEvent();

	Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int

	//Creating the histograms depending on MC (truth or reco) or Data
	createHistos();

	// TFile *outputFile = wk()->getOutputFile (_outputName);
	TFile *output = dynamic_cast<TFile *>(wk()->getOutputFile(_outputName));
     tree = new TTree ("tree", "tree");
     tree->SetDirectory (output);
//
     tree->Branch("XsectiontimesBR",   	        &xsectionbr);
//     tree->Branch("tot_weights", 	        &total_weights);
//     tree->Branch("analys_weight",           &weight_analyse);
//     tree->Branch("mass_yy", 	            &mass_yy);

	return StatusCode::SUCCESS;
}


EL::StatusCode HcAnalyse :: execute ()
{

    xAOD::TEvent* event = wk()->xaodEvent();
	int statSize=1;
	if(m_eventCounter!=0){
		double power=std::floor(log10(m_eventCounter));
		statSize=(int)std::pow(10.,power);
	}
	if(m_eventCounter%statSize==0) std::cout << "Event: " << m_eventCounter << std::endl;
	m_eventCounter++;


    //--------------------------------------------------------------------------------------------------
    //           Event information
    //--------------------------------------------------------------------------------------------------
    //EventInfo//////////////////////////////////
    const xAOD::EventInfo *eventInfo = 0;
    if (! event->retrieve(eventInfo, "EventInfo").isSuccess()){
    	Error("execute()", "Failed to retrieve event info collection. Exiting.");
    	return EL::StatusCode::FAILURE;}
    isMC = false; if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){  isMC = true;}
    //HGamEventInfo/////////////////////////////
    const xAOD::EventInfo *HgameventInfo = 0;
    if (! event->retrieve(HgameventInfo, "HGamEventInfo").isSuccess()) {
    	Error("execute()", "Failed to retrieve Hgam event info collection. Exiting.");
    	return EL::StatusCode::FAILURE;}
    const xAOD::EventInfo *HgamTrutheventInfo = 0;
  	if (!event->retrieve(HgamTrutheventInfo, "HGamTruthEventInfo").isSuccess()){
  	Error("execute()", "Failed to retrieve Hgam event info collection. Exiting.");
  	return EL::StatusCode::FAILURE;}


	//--------------------------------------------------------------------------------------------------
	//----------------------------------------------
	//            Initialize the containers
	//----------------------------------------------
	const xAOD::PhotonContainer* photons = 0;
	const xAOD::JetContainer* jets = 0;
	//const xAOD::JetContainer* PFlowjets = 0;
	const xAOD::ElectronContainer* electrons = 0;
	//const xAOD::MuonContainer* muons = 0;
 	const xAOD::TruthParticleContainer *truphotons = 0;

 	const xAOD::JetContainer *truthjets= 0;

    //Jets//Photons//Electrons//Muons
    EL_RETURN_CHECK("execute()",event->retrieve( jets, _reco_jet_collection.c_str() ));//pass the string of your favorite jet container
    //EL_RETURN_CHECK("execute()",event->retrieve( PFlowjets, "HGamAntiKt4EMPFlowJets_BTagging201903" ));//pass the string of your favorite jet container
    EL_RETURN_CHECK("execute()",event->retrieve( photons, "HGamPhotons" ));
    //--------------------------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------------------------
	// MAke truth histograms
//     if(isMC){
//   //   	TruthHGamEvent info///////////////////////
  //		const xAOD::EventInfo *HgamTrutheventInfo = 0;
  //		if (!event->retrieve(HgamTrutheventInfo, "HGamTruthEventInfo").isSuccess()){
//  		Error("execute()", "Failed to retrieve Hgam event info collection. Exiting.");
//  		return EL::StatusCode::FAILURE;}
//      	EL_RETURN_CHECK("execute()",event->retrieve( truphotons, "HGamTruthPhotons" )); if(_istruth)dotruth(truphotons,HgamTrutheventInfo);
// //

		// if(HgamTrutheventInfo->auxdataConst<char>("isFiducial")<1) {return EL::StatusCode::SUCCESS;};
//
// 		float mcweight;
// 		if(eventInfo->auxdataConst<float>("mcChannelNumber")==100000){mcweight= eventInfo->auxdataConst<float>("mcEventWeights");}else mcweight=1.0;
//
// 		HStore()->fillTH1F("N_j",HgamTrutheventInfo->auxdataConst<int>("N_j"),mcweight);
// 		HStore()->fillTH1F("N_j_ctag25",HgamTrutheventInfo->auxdataConst<int>("N_j_ctag25"),mcweight);
// 		HStore()->fillTH1F("dR_y_y",HgamTrutheventInfo->auxdataConst<float>("dR_y_y"),mcweight);
// 		HStore()->fillTH1F("pT_j1",HgamTrutheventInfo->auxdataConst<float>("pT_j1")* HC::invGeV,mcweight);
// 		HStore()->fillTH1F("pT_j2",HgamTrutheventInfo->auxdataConst<float>("pT_j2")* HC::invGeV,mcweight);
// 		HStore()->fillTH1F("pT_yy",HgamTrutheventInfo->auxdataConst<float>("pT_yy")* HC::invGeV,mcweight);
// 		HStore()->fillTH1F("y_yy",HgamTrutheventInfo->auxdataConst<float>("y_yy"),mcweight);
// 		HStore()->fillTH1F("pT_y1",HgamTrutheventInfo->auxdataConst<float>("pT_y1")* HC::invGeV,mcweight);
// 		HStore()->fillTH1F("pT_y2",HgamTrutheventInfo->auxdataConst<float>("pT_y2")* HC::invGeV,mcweight);

    //}
//     --------------------------------------------------------------------------------------------------

	// int isCjet = 0;
// 	int isCjet30 = 0;
// 	int isCjet35 = 0;
// 	int isCjet40 = 0;
// 	int isCjet50 = 0;

    //--------------------------------------------------------------------------------------------------
    //Do the cutflow of the NTuples  (pass the boolean - do_cutflow?)
    if(_do_Cutflow){
		fillCutFlow(CutEnum(ALL));
 		if(isMC){
     		//TruthHGamEvent info///////////////////////
 			const xAOD::EventInfo *HgamTrutheventInfo = 0;
 			if (!event->retrieve(HgamTrutheventInfo, "HGamTruthEventInfo").isSuccess()){
 			Error("execute()", "Failed to retrieve Hgam event info collection. Exiting.");
 			return EL::StatusCode::FAILURE;}
      		EL_RETURN_CHECK("execute()",event->retrieve( truthjets, "HGamAntiKt4TruthWZJets" ));

//
//       	for(auto truejet: *truthjets){
//       		    if(truejet->pt()*0.001 <25)continue;
//       		    if(fabs(truejet->eta())>2.5)continue;
//
//       		    //HStore()->fillTH1F("truejets_pt",truejet->pt()*0.001); HStore()->fillTH1F("truejets_eta",truejet->eta()*0.001);
//       		    //if(truejet->auxdataConst<int>("HadronConeExclTruthLabelID")==4){HStore()->fillTH1F("trueCjets_pt",truejet->pt()*0.001); HStore()->fillTH1F("trueCjets_eta",truejet->eta()*0.001);}
//
// 				if(truejet->auxdataConst<int>("HadronConeExclTruthLabelID")==4){
// 					isCjet++;
// 					if(truejet->pt()*0.001 >30)isCjet30++;
// 					if(truejet->pt()*0.001 >35)isCjet35++;
// 					if(truejet->pt()*0.001 >40)isCjet40++;
// 					if(truejet->pt()*0.001 >50)isCjet50++;
// 				}
//       		}
      	}

		//if(HgameventInfo->auxdataConst<char>("isPassedPreselection")<1) {return EL::StatusCode::SUCCESS;}fillCutFlow(CutEnum(PRESEL));
		//2. passed selection
		//if(HgameventInfo->auxdataConst<char>("isPassed")!=1) {return EL::StatusCode::SUCCESS;}fillCutFlow(CutEnum(SEL));

    }//end of cutflow
    //--------------------------------------------------------------------------------------------------


    //Checking the weights which include the InitialWeights*SF
    WeightsScales(HgameventInfo);


  //   auto leading_photon = (*photons)[0];
//     auto subleading_photon = (*photons)[1];
//     if(_ispTcut){ //we cut on the pt cut for the photon which is 25 GeV and 22
//        if((leading_photon->pt() < _phpTcut * HC::GeV ) || (subleading_photon->pt() < _phpTcut * HC::GeV))return EL::StatusCode::SUCCESS;
//      }
//      //Cutflow
//      fillCutFlow(CutEnum(PHPT));
//


//      HStore()->fillTH1F("N_j_all",n_jets,weight_analyse);
//
//
//     //Now try to do the selection with Nj_30
//     int nj30 = 0;
//     for(auto recoj: *jets){
//     	if(recoj->pt() < 30 * HC::GeV )continue;
//     	nj30++;
//     }
//
//     HStore()->fillTH1F("N_j30",nj30,weight_analyse30);


	 //Here is the mass before the jet cut
	mass_gev = HgameventInfo->auxdataConst<float>("m_yy") * HC::invGeV;
	//mass_yy_pv = HgameventInfo->auxdataConst<float>("m_yy_hardestVertex") * HC::invGeV;
	//mass_yy_true = HgameventInfo->auxdataConst<float>("m_yy_truthVertex") * HC::invGeV;

//
// 	HStore()->fillTH1F("m_yy_nojet",mass_gev,total_weights);
//     if(mass_gev>130.)return EL::StatusCode::SUCCESS;
//     if(mass_gev<120.)return EL::StatusCode::SUCCESS;
//     HStore()->fillTH1F("m_yy_nojet_inside",mass_gev,total_weights);

    xsectionbr = HgameventInfo->auxdataConst<float>("crossSectionBRfilterEff");

   // if(HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid")>0)return EL::StatusCode::SUCCESS;

  //   int ismatched;
//     if(mass_gev == mass_yy_pv){ismatched = 1;}
//     else {ismatched = 0;}
//
//     HStore()->fillTH1F("NN_vs_PV",ismatched);
//     HStore()->fillTH1F("m_yy_pv",mass_yy_pv);
//     HStore()->fillTH1F("m_yy_nn",mass_gev);
//     HStore()->fillTH1F("m_yy_true",mass_yy_true);
//
//
//     int match_nn_true;
//     int match_pv_true;
//
//     if(mass_gev == mass_yy_true){match_nn_true = 1;}
//     else {match_nn_true=0;}
//
//     if(mass_yy_pv == mass_yy_true){match_pv_true = 1;}
//     else {match_pv_true=0;}


//    HStore()->fillTH1F("NN_vs_TRUE",match_nn_true);
//    HStore()->fillTH1F("PV_vs_TRUE",match_pv_true);

// 	if(jets->size()>0){
//     	auto leadjet = (*jets)[0];
// 		if(fabs(leadjet->eta())<2.5){HStore()->fillTH1F("customlead_jet_pt",leadjet->pt()*0.001);}
// 	}
//
// 	if(PFlowjets->size()>0){
// 		auto pflowleadjet = (*PFlowjets)[0];
// 		if(fabs(pflowleadjet->eta())<2.5){HStore()->fillTH1F("pflowlead_jet_pt",pflowleadjet->pt()*0.001);}
// 	}

    //Cutting on number og jets in the events
  //   int n_jets = jets->size();
//     if(n_jets < 1 )return EL::StatusCode::SUCCESS;
//     fillCutFlow(CutEnum(JET1));
//     //--------------------------------------------------------------------------------------------------

    //the masses calculated by hand and m_yy form Ntup are identical! -------------------------(checked)
    //Fill the histo (without and with jvt weight) for the at least 1 jet
   // HStore()->fillTH1F("m_yy_1j_jvt",mass_gev,weight_analyse);

    HStore()->fillTH1F("N_cjets_weighted",HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid"),HgameventInfo->auxdataConst<float>("weightInitial"));
    HStore()->fillTH1F("N_cjets",HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid"));
    //--------------------------------------------------------------------------------------------------

//
//
//     //Cehcing form the Ntuple the m_yy full with weights
//     if(HgameventInfo->auxdataConst<int>("Hc_Atleast1jisloose")==0){
//         HStore()->fillTH1F("m_yy_full_nonctag",mass_gev,weight_analyse_jvt_ctagnom);
// //        HStore()->fillTH1F("m_yy_full_nonctag_jvt",mass_gev,weight_analyse_jvt);
// //	    HStore()->fillTH1F("m_yy_full_nonctag_light",mass_gev,weight_analyse_jvt_light);
//
//
//      	if(HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid")<1){
//      		HStore()->fillTH1F("m_yy_full_nonctag_bkg",mass_gev,weight_analyse_jvt_ctagnom);
// // 			HStore()->fillTH1F("m_yy_full_nonctag_bkg_jvt",mass_gev,weight_analyse_jvt);
// // 			HStore()->fillTH1F("m_yy_full_nonctag_bkg_light",mass_gev,weight_analyse_jvt_light);
//  		}
// //
//      	if(HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid")>0){
//      		HStore()->fillTH1F("m_yy_full_nonctag_sig",mass_gev,weight_analyse_jvt_ctagnom);
// // 			HStore()->fillTH1F("m_yy_full_nonctag_sig_jvt",mass_gev,weight_analyse_jvt);
// // 			HStore()->fillTH1F("m_yy_full_nonctag_sig_light",mass_gev,weight_analyse_jvt_light);
//      	}


//      	if(isCjet==1){HStore()->fillTH1F("m_yy_full_nonctag_sig_1j",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet>1){HStore()->fillTH1F("m_yy_full_nonctag_sig_1jplus",mass_gev,weight_analyse_jvt_ctagnom);}
//
// 		if(isCjet30>0){HStore()->fillTH1F("m_yy_full_nonctag_sig30",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet30<1) {HStore()->fillTH1F("m_yy_full_nonctag_bkg30",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet35>0){HStore()->fillTH1F("m_yy_full_nonctag_sig35",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet35<1) {HStore()->fillTH1F("m_yy_full_nonctag_bkg35",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet40>0){HStore()->fillTH1F("m_yy_full_nonctag_sig40",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet40<1) {HStore()->fillTH1F("m_yy_full_nonctag_bkg40",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet50>0){HStore()->fillTH1F("m_yy_full_nonctag_sig50",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet50<1) {HStore()->fillTH1F("m_yy_full_nonctag_bkg50",mass_gev,weight_analyse_jvt_ctagnom);}

 //   }

//
//      if(HgameventInfo->auxdataConst<int>("Hc_Atleast1jisloose")==1){
//         HStore()->fillTH1F("m_yy_full_ctag",mass_gev,weight_analyse_jvt_ctagnom);
// //		HStore()->fillTH1F("m_yy_full_ctag_jvt",mass_gev,weight_analyse_jvt);
// //		HStore()->fillTH1F("m_yy_full_ctag_light",mass_gev,weight_analyse_jvt_light);
//
//  		if(HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid")<1){
//  	    	HStore()->fillTH1F("m_yy_full_ctag_bkg",mass_gev,weight_analyse_jvt_ctagnom);
// // 			HStore()->fillTH1F("m_yy_full_ctag_bkg_jvt",mass_gev,weight_analyse_jvt);
// // 			HStore()->fillTH1F("m_yy_full_ctag_bkg_light",mass_gev,weight_analyse_jvt_light);
//      	}
//      	if(HgamTrutheventInfo->auxdataConst<int>("N_j_cjet_had_fid")>0){
//      		HStore()->fillTH1F("m_yy_full_ctag_sig",mass_gev,weight_analyse_jvt_ctagnom);
// // 			HStore()->fillTH1F("m_yy_full_ctag_sig_jvt",mass_gev,weight_analyse_jvt);
// // 			HStore()->fillTH1F("m_yy_full_ctag_sig_light",mass_gev,weight_analyse_jvt_light);
//      	}

     // 	if(isCjet==1){HStore()->fillTH1F("m_yy_full_ctag_sig_1j",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet>1){HStore()->fillTH1F("m_yy_full_ctag_sig_1jplus",mass_gev,weight_analyse_jvt_ctagnom);}
//
//
// 		if(isCjet30>0){HStore()->fillTH1F("m_yy_full_ctag_sig30",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet30<1) {HStore()->fillTH1F("m_yy_full_ctag_bkg30",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet35>0){HStore()->fillTH1F("m_yy_full_ctag_sig35",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet35<1) {HStore()->fillTH1F("m_yy_full_ctag_bkg35",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet40>0){HStore()->fillTH1F("m_yy_full_ctag_sig40",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet40<1) {HStore()->fillTH1F("m_yy_full_ctag_bkg40",mass_gev,weight_analyse_jvt_ctagnom);}
// 		if(isCjet50>0){HStore()->fillTH1F("m_yy_full_ctag_sig50",mass_gev,weight_analyse_jvt_ctagnom);}else if(isCjet50<1) {HStore()->fillTH1F("m_yy_full_ctag_bkg50",mass_gev,weight_analyse_jvt_ctagnom);}

//   }



	// int ismatched_nonctag;
// 	int ismatched_ctagged;
//
// 	int ismatched_nonctag_nn_true;
// 	int ismatched_ctagged_nn_true;
//
// 	int ismatched_nonctag_pv_true;
// 	int ismatched_ctagged_pv_true;
//
// 	if(HgameventInfo->auxdataConst<int>("Hc_Atleast1jisloose")==0){
// 		if((mass_gev == mass_yy_pv) ) {ismatched_nonctag = 1;}
// 		else {ismatched_nonctag = 0;}
//
// 		if((mass_gev == mass_yy_true) ) {ismatched_nonctag_nn_true= 1;}
// 		else {ismatched_nonctag_nn_true = 0;}
//
// 		if((mass_yy_pv == mass_yy_true) ) {ismatched_nonctag_pv_true = 1;}
// 		else {ismatched_nonctag_pv_true = 0;}
// 	}
//
// 	if(HgameventInfo->auxdataConst<int>("Hc_Atleast1jisloose")==1){
// 		if((mass_gev == mass_yy_pv) ) {ismatched_ctagged = 1;}
// 		else {ismatched_ctagged = 0;}
//
// 		if((mass_gev == mass_yy_true) ) {ismatched_ctagged_nn_true= 1;}
// 		else {ismatched_ctagged_nn_true = 0;}
//
// 		if((mass_yy_pv == mass_yy_true) ) {ismatched_ctagged_pv_true = 1;}
// 		else {ismatched_ctagged_pv_true = 0;}
// 	}
//
// 	HStore()->fillTH1F("NN_vs_PV_ctagged",ismatched_ctagged);
// 	HStore()->fillTH1F("NN_vs_PV_nonctag",ismatched_nonctag);
//
// 	HStore()->fillTH1F("NN_vs_TRUE_ctagged",ismatched_ctagged_nn_true);
// 	HStore()->fillTH1F("NN_vs_TRUE_nonctag",ismatched_nonctag_nn_true);
//
// 	HStore()->fillTH1F("PV_vs_TRUE_ctagged",ismatched_ctagged_pv_true);
// 	HStore()->fillTH1F("PV_vs_TRUE_nonctag",ismatched_nonctag_pv_true);

// 	if(HgameventInfo->auxdataConst<int>("Hcgam_Atleast1jisloose")==1 && HgameventInfo->auxdataConst<int>("Hcgam_CountCtruthFromCharm")==1){
// 		 HStore()->fillTH1F("m_yy_ctag_jvt",mass_gev,weight_analyse);
// 	}
//
// 	if(HgameventInfo->auxdataConst<int>("Hcgam_Atleast1jisloose")==0 && HgameventInfo->auxdataConst<int>("Hcgam_CountCtruthFromCharm")==1){
// 		 HStore()->fillTH1F("m_yy_nonctag_jvt",mass_gev,weight_analyse);
// 	}
    //Check the jsets with different pT and /or check the optimization
//    checkJets(HgameventInfo, jets, photons, weight_analyse_jvt, weight_analyse_jvt_ctagnom, weight_analyse_jvt_light);


    //Checking he distributions of photons
    //if(_isreco)doreco(photons,jets);


     tree->Fill();
  return StatusCode::SUCCESS;
}



EL::StatusCode HcAnalyse :: postExecute (){ return EL::StatusCode::SUCCESS;}


EL::StatusCode HcAnalyse :: finalize ()
{
  if(_do_Cutflow){
  	printf("\nEvent selection cut flow:\n");
  	hist_out=getCutFlowHisto();
  	printCutFlowHisto(hist_out,0);
  }

  SafeDelete(m_hStore);
  return StatusCode::SUCCESS;
}

EL::StatusCode HcAnalyse :: histFinalize ()
{
  std::cout<<"Events = "<<m_eventCounter<<std::endl;
  return EL::StatusCode::SUCCESS;
}

//Some functions used in the processing
//--------------------------------------------------------------------------------------------------
void HcAnalyse::printCutFlowHisto(TH1F *h, int Ndecimals)
{
  TString format("  %-24s%10.");
  format += Ndecimals;
  format += "f%11.2f%%%11.2f%%\n";
  int all_bin = h->FindBin(0);
  printf("  %-24s%10s%12s%12s\n", "Event selection", "Nevents", "Cut rej.", "Tot. eff.");

  for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
    double ntot = h->GetBinContent(all_bin), n = h->GetBinContent(bin), nprev = h->GetBinContent(bin - 1);
    TString cutName(h->GetXaxis()->GetBinLabel(bin));
    cutName.ReplaceAll("#it{m}_{#gamma#gamma}", "m_yy");

    if (bin == 1 || nprev == 0 || n == nprev)
    { printf(format.Data(), cutName.Data(), n, -1e-10, n / ntot * 100); }
    else // if the cut does something, print more information
    { printf(format.Data(), cutName.Data(), n, (n - nprev) / nprev * 100, n / ntot * 100); }
  }
}

//--------------------------------------------------------------------------------------------------
//Dealing with the weights:
void HcAnalyse::WeightsScales(const xAOD::EventInfo* HgameventInfo)
{

    //==1. - initial weight - to be used for truth histos (mcWeight * pileupWeight * vertexWeight)-----(checked)
    initial_weight = HgameventInfo->auxdataConst<float>("weightInitial");

    float final_weight  = HgameventInfo->auxdataConst<float>("weight");
     //float jvt25 = HgameventInfo->auxdataConst<float>("weightJvt");//weight for the jets above 25 GeV
    // float jvt30 = HgameventInfo->auxdataConst<float>("weightJvt_30");//weight for the jets above 25 GeV
    //float jvt50 = HgameventInfo->auxdataConst<float>("weightJvt_50");//weight for the jets above 25 GeV
    //float ctagSF    = HgameventInfo->auxdataConst<float>("Hc_weight_ctag_e");//weight for the jets above 25 GeV
    //float nonctagSF = HgameventInfo->auxdataConst<float>("Hc_weight_nonctag_e");//weight for the jets above 25 GeV

    //weight_sf - already includes the weights from Reco*Trig for muons--------------------(checked)
    float weight_sf = HgameventInfo->auxdataConst<float>("weightSF");
    //==2. - total weight - to be used for the reco histos - the same as final weight
    total_weights = initial_weight * weight_sf; //tot_weight_nosf = initial_weight * jvt25;
    //Fill the histos
    HStore()->fillTH1F("total_weights_calc",total_weights);
    HStore()->fillTH1F("total_weights_Ntup",final_weight);//the same -->checked

    //float jvt  = HgameventInfo->auxdataConst<float>("Hc_weightjvt");
    //float ctag = HgameventInfo->auxdataConst<float>("Hc_weightCtag");
    float light= 1; //HgameventInfo->auxdataConst<float>("Hc_weightCtagLight");

    weight_analyse_jvt         = total_weights;// * jvt25;
   //weight_analyse_jvt         = total_weights * jvt;
    weight_analyse_jvt_ctagnom = total_weights; //* jvt * ctag;
    weight_analyse_jvt_light   = total_weights; //* jvt * light;


    //With the JVT weights
    //==3. - Calculate the total weight with jvt
   // weight_analyse = total_weights * jvt25;
    //weight_analyse30 = total_weights * jvt30;
    //weight_analyse50 = total_weights; //* jvt50;

   // weight_ctag = weight_analyse*ctagSF;
   // weight_nonctag = weight_analyse*nonctagSF;


    //Fill the histo
    HStore()->fillTH1F("jvt_plus_total",weight_analyse_jvt_ctagnom);
    HStore()->fillTH1F("jvt_weight",weight_analyse_jvt);



}
//--------------------------------------------------------------------------------------------------
//Dealing with jet and ctag optimization
EL::StatusCode HcAnalyse::checkJets(const xAOD::EventInfo* HgameventInfo,const xAOD::JetContainer *jets, const xAOD::PhotonContainer *photon, float weightjvt, float weightjvt_ctag, float weightjvt_light)
{
//
		 auto leadjet = (*jets)[0];
		//Lets make the photon cuts
		auto pho_lead = (*photon)[0];
	    auto pho_sublead = (*photon)[1];
		TLorentzVector pho1 = pho_lead->p4(), pho2 = pho_sublead->p4();
	    pho1 *= HC::invGeV; pho2 *= HC::invGeV; TLorentzVector Higgs = pho1+pho2;

		//std::cout<<"leading ptohont pt "<<pho_lead->pt()/1000.<<std::endl;
		//std::cout<<"subleading ptohont pt "<<pho_sublead->pt()/1000.<<std::endl;
	    TLorentzVector jet = leadjet->p4();
	    jet *= HC::invGeV;
	    TLorentzVector jgamgam = jet + Higgs;
//
// 		int countcjets_b70 = 0;
// 		int countcjets_nob_loose  = 0;
// 		int countcjets_nob_looser = 0;
// 		int countcjets_nob_ofic   = 0;
// 		//int countcjets_nob_tight = 0;
//
// 		//int countcjets_b70_j30 = 0;
//
// 		int countcjets_nob_loose_j30  = 0;
// 		int countcjets_nob_looser_j30 = 0;
// 		int countcjets_nob_ofic_j30   = 0;
//
//
// 		//int countcjets_nob_tight_j30 = 0;
//
//
// 		//int counterLoose = -1, counterTight = -1;
// 	    //int firstCtagLoose = -3, firstCtagTight = -3;
// 	    int drcounts = -1;
// 	    int firstdrcounts = -2;
//
//
//
// 	    int counterLoose = -1, counterLooser = -1, counterOfic = -1;
// 	    int firstCtagLoose = -3, firstCtagLooser = -3, firstCtagOfic = -3;
// 	    int passing2 = 0;
//
// 		int counts=0;
// 	    for(auto recoj: *jets){
// 	    	if(fabs(recoj->eta())>2.5)continue;
// 	    	if(recoj->pt()< 25 * HC::GeV )continue;
// 	    	//HStore()->fillTH1F("alljet_pt_25",recoj->pt() * HC::invGeV, weightjvt);
//
// 	    	counts++;
// 	    	if(counts>=3)break;
// 			//for the selection
// 			if(firstCtagLoose  ==-3) {firstCtagLoose  = -2;}
// 			if(firstCtagLooser ==-3) {firstCtagLooser = -2;}
// 			if(firstCtagOfic   ==-3) {firstCtagOfic   = -2;}
// 			//if(firstCtagTight ==-3) {firstCtagTight = -2;}
// 			if(firstdrcounts ==-2){firstdrcounts = -1;}
// 			bool isCjet = false, isBjet = false, isLjet = false;
//     		bool mv2_70j = false;
//
//     		//DL1
// 			//double pu_j = recoj->auxdataConst<double>("DL1_pu");
// 			//double pc_j = recoj->auxdataConst<double>("DL1_pc");
// 			//double pb_j = recoj->auxdataConst<double>("DL1_pb");
//
// 			double pu_jr = recoj->auxdataConst<double>("DL1r_pu");
// 			double pc_jr = recoj->auxdataConst<double>("DL1r_pc");
// 			double pb_jr = recoj->auxdataConst<double>("DL1r_pb");
//
// 			counterLoose++;
// 			counterLooser++;
// 			counterOfic++;
// 			//counterTight++;
//
//
// 			drcounts++;
//
// 			//Dl
// 			//double dl1_discriminant_DL1lose   = discriminant(pc_j,pb_j,pu_j,0.0225);
// 			double dl1_discriminant_DL1loser  = discriminant(pc_jr,pb_jr,pu_jr,0.0235);
// 			double dl1_discriminant_DL1Ofic   = discriminant(pc_jr,pb_jr,pu_jr,0.2);
// 			//double dl1_discriminant_DL1tight = discriminant(pc_j,pb_j,pu_j,0.08);
//
//
// 			if(HC::DR(recoj,pho_lead) < _dR_cut || HC::DR(recoj,pho_sublead) < _dR_cut) continue;
//
// 			if(firstCtagLoose  ==-2) {firstCtagLoose  = -1;}
// 			if(firstCtagLooser ==-2) {firstCtagLooser = -1;}
// 			if(firstCtagOfic   ==-2) {firstCtagOfic   = -1;}
// 			//if(firstdrcounts   ==-2)  {firstdrcounts  = -1;}
// 			//if(firstCtagTight ==-2) {firstCtagTight = -1;}
//
// 			firstdrcounts = drcounts;
//
// 			//if(drcounts ==0){
// 			//	HStore()->fillTH1F("dR_j1y1",HC::DR(recoj,pho_lead));
// 			//	HStore()->fillTH1F("dR_j1y2",HC::DR(recoj,pho_sublead));
//
// 			//}
// 			//if(firstdrcounts ==-1){firstdrcounts =drcounts;}
//
// 			if(dl1_discriminant_DL1Ofic < 1  && drcounts==0){HStore()->fillTH1F("TruthLabel_NonCtagTight_lead",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));}
// 			if(dl1_discriminant_DL1Ofic < 1  && drcounts==1){HStore()->fillTH1F("TruthLabel_NonCtagTight_sublead",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));}

			//if(mv2_70j==false && dl1_discriminant_DL1tight > 1.3){countcjets_b70 ++;}
			//if(dl1_discriminant_DL1lose > 1.44){
			//	countcjets_nob_loose ++;
			//	//HStore()->fillTH1F("TruthLabel_AllCtagLoose_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
			//}

			//if(dl1_discriminant_DL1loser > 1.44){
			//	countcjets_nob_looser ++;
			//	//HStore()->fillTH1F("TruthLabel_AllCtagLoose_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
			//}

			//if(dl1_discriminant_DL1Ofic > 1.){
			//	countcjets_nob_ofic ++;
			//	//HStore()->fillTH1F("TruthLabel_AllCtagLoose_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
			//}

			//if(dl1_discriminant_DL1tight > 1.3){
			//	countcjets_nob_tight ++;
			//	//HStore()->fillTH1F("TruthLabel_AllCtagTight_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
			//}

			//if(dl1_discriminant_DL1lose > 1.44 && firstCtagLoose==-1){
			//	firstCtagLoose = counterLoose;
				//HStore()->fillTH1F("TruthLabel_FirstCtagLoose_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
			//}

			// if(dl1_discriminant_DL1loser > 1.53 && firstCtagLooser==-1){
// 				firstCtagLooser = counterLooser;
// 				//HStore()->fillTH1F("TruthLabel_FirstCtagTight_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
// 			}

	// 		if(dl1_discriminant_DL1Ofic > 1 && firstCtagOfic==-1){
// 				firstCtagOfic = counterOfic;
// 				HStore()->fillTH1F("TruthLabel_FirstCtagTight_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
// 			}
//
// 		}

		// if(firstCtagOfic==-1){
// 		if (jets->size()>0){
// 		HStore()->fillTH1F("TruthLabel_NonCtagTight_lead",(*jets)[0]->auxdataConst<int>("HadronConeExclTruthLabelID"));
// 		HStore()->fillTH1F("TruthLabel_NonCtagTight_sublead",(*jets)[1]->auxdataConst<int>("HadronConeExclTruthLabelID"));
// 		}
// 		}

		//HStore()->fillTH1F("Jet_index_drcounts",firstdrcounts);



		//HStore()->fillTH1F("Jet_index_1ctagLoose_j25",firstCtagLoose);
		//HStore()->fillTH1F("Jet_index_1ctagLooser_j25",firstCtagLooser);
		//HStore()->fillTH1F("Jet_index_1ctagOfic_j25",firstCtagOfic);


		//HStore()->fillTH1F("Jet_index_1ctagTight_j25",firstCtagTight);

		//if(firstCtagLoose>=-1){HStore()->fillTH1F("m_yy_allj_25gev",mass_gev,weightjvt);}


		//if(firstCtagLoose==-1){HStore()->fillTH1F("m_yy_allj_25gev",mass_gev,weightjvt);}

		//if( firstCtagLoose>=0  && firstCtagLoose<2 ) {HStore()->fillTH1F("m_yy_loose",mass_gev,weightjvt);}  else if(firstCtagLoose==-1)  {HStore()->fillTH1F("m_yy_NOloose",mass_gev,weightjvt);}  //&& countcjets_nob_loose>0 ) {HStore()->fillTH1F("m_yy_loose",mass_gev,weightjvt);} else {HStore()->fillTH1F("m_yy_NOloose",mass_gev,weightjvt);}
		//if( firstCtagLooser>=0 && firstCtagLooser<2) {HStore()->fillTH1F("m_yy_looser",mass_gev,weightjvt);} else if(firstCtagLooser==-1) {HStore()->fillTH1F("m_yy_NOlooser",mass_gev,weightjvt);} //&& countcjets_nob_looser>0 )

//
// 		//if( firstCtagOfic>=0   && firstCtagOfic<2)   {HStore()->fillTH1F("m_yy_ofic_jvt",mass_gev,weightjvt);}   else if(firstCtagOfic==-1)   {HStore()->fillTH1F("m_yy_NOofic_jvt",mass_gev,weightjvt);}//&& countcjets_nob_ofic>0 )
// //
//  		if( firstCtagOfic>=0   && firstCtagOfic<2)   {HStore()->fillTH1F("m_yy_ofic_jvt_ctag",mass_gev,weightjvt_ctag);}   else if(firstCtagOfic==-1) {HStore()->fillTH1F("m_yy_NOofic_jvt_ctag",mass_gev,weightjvt_ctag);}//&& countcjets_nob_ofic>0 )
// //
// // 		if( firstCtagOfic>=0   && firstCtagOfic<2)   {HStore()->fillTH1F("m_yy_ofic_jvt_light",mass_gev,weightjvt_light);}   else if(firstCtagOfic==-1)   {HStore()->fillTH1F("m_yy_NOofic_jvt_light",mass_gev,weightjvt_light);}//&& countcjets_nob_ofic>0 )
//
//        	if(passing2 ==1){HStore()->fillTH1F("m_yy_ofic_2taggedCtag",mass_gev,weightjvt);} else if( firstCtagOfic==0){HStore()->fillTH1F("m_yy_ofic_leadCtag",mass_gev,weightjvt);} else if( firstCtagOfic==1){HStore()->fillTH1F("m_yy_ofic_subleadCtag",mass_gev,weightjvt);}
//

		//if( (firstCtagTight>=0 && firstCtagTight<2) && countcjets_nob_tight>0){HStore()->fillTH1F("m_yy_tight",mass_gev,weightjvt);}
		//if(countcjets_b70>0){HStore()->fillTH1F("m_yy_tight_nob",mass_gev,weightjvt);}


		//if(firstdrcounts>=0 && firstdrcounts<2){HStore()->fillTH1F("m_yy_allj_25gev",mass_gev,weightjvt);}




// 		int n_btagged_jets = 0; // initialize the counter for b-tagged jets
// 		int n_nonbtagged_jets = 0; // initialize the counter for non-btagged jets
// 		int jetcount = 0;

        int n_btagged_jets; // initialize the counter for b-tagged jets
		int n_nonbtagged_jets; // initialize the counter for non-btagged jets


//
// 		for(auto recoj: *jets){
//     		if(fabs(recoj->eta()) > 2.5) continue;
//     		if(recoj->pt() < 25 * HC::GeV) continue;
//     		jetcount++;
//
//     		if(jetcount>=3)break;
//
//     		if(HC::DR(recoj,pho_lead) < _dR_cut || HC::DR(recoj,pho_sublead) < _dR_cut) continue;
//
//
//     		double pu_jr = recoj->auxdataConst<double>("DL1r_pu");
//     		double pc_jr = recoj->auxdataConst<double>("DL1r_pc");
//     		double pb_jr = recoj->auxdataConst<double>("DL1r_pb");
//
//     		double dl1_discriminant_DL1Ofic = discriminant(pc_jr, pb_jr, pu_jr, 0.2);
//
//     		if(dl1_discriminant_DL1Ofic > 1) {
//         		n_btagged_jets++;
//         		HStore()->fillTH1F("TruthLabel_CtagTight_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
//
//     		}
//     		else if(dl1_discriminant_DL1Ofic < 1 ) {
//         		n_nonbtagged_jets++;
//         		HStore()->fillTH1F("TruthLabel_NonCtagTight_j25",recoj->auxdataConst<int>("HadronConeExclTruthLabelID"));
//     		}
//
//
//
//         }

//     	if(n_btagged_jets > 0 ) { HStore()->fillTH1F("m_yy_ofic_jvt_ctag",mass_gev,weightjvt_ctag);}
//     	else if(n_nonbtagged_jets>0) {HStore()->fillTH1F("m_yy_NOofic_jvt_ctag",mass_gev,weightjvt_ctag);}
//
//
//     	HStore()->fillTH1F("N_NoNctagged",n_nonbtagged_jets);
//     	HStore()->fillTH1F("N_Ctagged",n_btagged_jets);



    	if(_Dooptimization){
    		for(int i = 0; i < bins; i++){//dl cut loop loop
    			float score = count[i];
				for(int k = 0; k < fbins; k++){//fraction loop
					float bfrac = fcount[k];
					n_btagged_jets = 0; // initialize the counter for b-tagged jets
		            n_nonbtagged_jets = 0; // initialize the counter for non-btagged jets
		            int jetcount = 0;

    				for(auto recoj: *jets){
    					if(fabs(recoj->eta()) > 2.5) continue;
    					if(recoj->pt() < 25 * HC::GeV) continue;
    					jetcount++;
    					if(jetcount>=3)break;

    					if(HC::DR(recoj,pho_lead) < _dR_cut || HC::DR(recoj,pho_sublead) < _dR_cut) continue;


    					double pu_jr = recoj->auxdataConst<double>("DL1r_pu");
    					double pc_jr = recoj->auxdataConst<double>("DL1r_pc");
    					double pb_jr = recoj->auxdataConst<double>("DL1r_pb");

    					double dl1_discriminant_DL1Ofic = discriminant(pc_jr, pb_jr, pu_jr,bfrac);

    					if(dl1_discriminant_DL1Ofic > score){
    						n_btagged_jets++;
        				}
    			        else if(dl1_discriminant_DL1Ofic < score ){
        					n_nonbtagged_jets++;
    					}

					}//end of jets

				    if(n_btagged_jets > 0 ) {masses_mv60[i][k]->Fill(mass_gev,weight_analyse_jvt);}
    				else if(n_nonbtagged_jets > 0) {masses_mv60nonc[i][k]->Fill(mass_gev,weight_analyse_jvt);}

        		}//end of frac
        	}//end of dl1
		}//end of if statenment


return EL::StatusCode::SUCCESS;

}

//--------------------------------------------------------------------------------------------------
////////////////////////////////////////////
//Cutflow hist maker Inherited from the HGamTool////////////////
////////////////////////////////////////////
TH1F *HcAnalyse::makeCutFlowHisto(int id, TString suffix)
{
   // const char ids = "something";
   // TString name = Form("CutFlow_%s",ids);
    int Ncuts = s_cutDescs.size();
    TH1F *h = new TH1F(suffix, suffix, Ncuts, 0, Ncuts);
    for (int bin = 1; bin <= Ncuts; ++bin)
    { h->GetXaxis()->SetBinLabel(bin, s_cutDescs[bin - 1]); }
    wk()->addOutput(h);
	return h;
}

void HcAnalyse::fillCutFlow(CutEnum cut)
{
    getCutFlowHisto()->Fill(cut);
}


//--------------------------------------------------------------------------------------------------
//Discriminant function//very easy
double HcAnalyse::discriminant(double pc, double pb, double pu, double fraction)
{
    double ctag_discriminant = TMath::Log( pc / ( (pb * fraction) + (1 - fraction)*pu ) );
    return ctag_discriminant;
}


//--------------------------------------------------------------------------------------------------
//creating the histos
EL::StatusCode HcAnalyse::createHistos()
{

		//fiducial
		HStore()->createTH1F("N_j",8,-0.5,7.5,";#it{N}_{j}");
		HStore()->createTH1F("N_ctag25",8,-0.5,7.5,";#it{N}_{c}");
		HStore()->createTH1F("dR_y_y",100,0.,10.,";#it{dR}^{#gamma#gamma}");
		HStore()->createTH1F("y_yy",20, -2.5, 2.5,";#it{y}^{#gamma#gamma}");
		HStore()->createTH1F("pT_y1",100, 0., 200.,";#it{pT}^{#gamma1} [GeV]");
		HStore()->createTH1F("pT_y2",100, 0., 200.,";#it{pT}^{#gamma2} [GeV]");
		HStore()->createTH1F("pT_j1",100, 0., 200.,";#it{pT}^{j1} [GeV]");
		HStore()->createTH1F("pT_j2",100, 0., 200.,";#it{pT}^{j2} [GeV]");
		HStore()->createTH1F("pT_yy",100, 0., 200.,";#it{pT}^{#gamma#gamma} [GeV]");



	if(_isreco){
		//For the plots of reconstructed photons
		HStore()->createTH1F("myy", 60, 100., 160.,";#it{m}_{#gamma#gamma} [GeV]");
		HStore()->createTH1F("ptyy",200, 0., 200.,";#it{p}^{#gamma#gamma}_{T} [GeV]");
		HStore()->createTH1F("yyy", 25, -2.5, 2.5,";#it{y}^{#gamma#gamma}");
		HStore()->createTH1F("dryy", 100, 0., 10.,";#it{dr}^{#gamma#gamma}");
		HStore()->createTH1F("dyyy", 50, -5., 5.,";#it{dy}^{#gamma#gamma}");
		HStore()->createTH1F("dphiyy", 64, -8., 8.,";#it{dphi}^{#gamma#gamma}");
		//singles
		HStore()->createTH1F("y1pt",200, 0., 200.,";#it{p}^{#gamma}_{T} [GeV]");
		HStore()->createTH1F("y2pt",200, 0., 200.,";#it{p}^{#gamma}_{T} [GeV]");
		HStore()->createTH1F("y1eta",50,-2.5, 2.5,";#it{#eta}^{#gamma}");
		HStore()->createTH1F("y2eta",50,-2.5, 2.5,";#it{#eta}^{#gamma}");
		HStore()->createTH1F("y1phi",30,-3.2, 3.2,";#it{#phi}^{#gamma}");
		HStore()->createTH1F("y2phi",50,-3.2, 3.2,";#it{#phi}^{#gamma}");

		HStore()->createTH1F("dRjgam1", 100, 0., 10.,";#it{dr}^{#gamma{1}j}");
		HStore()->createTH1F("dRjgam2", 100, 0., 10.,";#it{dr}^{#gammф{2}j}");
		HStore()->createTH1F("dphijgam1", 64, -8., 8.,";#it{dphi}^{#gamma_{1}j}");
		HStore()->createTH1F("dphijgam2", 64, -8., 8.,";#it{dphi}^{#gamma_{2}j}");

		HStore()->createTH1F("jpt",400,0.,400.,";#it{p^{j}_{T}");
	    HStore()->createTH1F("jy",25, -2.5, 2.5,";#it{y}^{j}");
	    HStore()->createTH1F("jeta",25, -2.5, 2.5,";#it{#eta}^{j}");

		HStore()->createTH1F("jgamgam1_pt",400,0.,400.,";#it{p^{j#gamma#gamma}_{T}");
		HStore()->createTH1F("jgamgam1_m",400,0.,400.,";#it{M^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_phi",50,-3.2, 3.2,";#it{#phi^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_eta",25, -2.5, 2.5,";#it{#eta^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_y",25, -2.5, 2.5,";#it{y^{j#gamma#gamma}");

		HStore()->createTH1F("jh_dR", 100, 0., 10.,";#it{dr}^{#gamma{1}j}");
		HStore()->createTH1F("jh_dphi",64, -8., 8.,";#it{dr}^{#gamma{1}j}");
		HStore()->createTH1F("jh_dy", 50, -5., 5.,";#it{dr}^{#gamma{1}j}");

		HStore()->createTH1F("dRjgam1_30", 100, 0., 10.,";#it{dr}^{#gamma{1}j}");
		HStore()->createTH1F("dRjgam2_30", 100, 0., 10.,";#it{dr}^{#gammф{2}j}");
		HStore()->createTH1F("dphijgam1_30", 64, -8., 8.,";#it{dphi}^{#gamma_{1}j}");
		HStore()->createTH1F("dphijgam2_30", 64, -8., 8.,";#it{dphi}^{#gamma_{2}j}");



	    HStore()->createTH1F("jeta_30",25, -2.5, 2.5,";#it{#eta}^{j}");


		HStore()->createTH1F("jgamgam1_m_30",400,0.,400.,";#it{M^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_pt_30",400,0.,400.,";#it{M^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_phi_30",50,-3.2, 3.2,";#it{#phi^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_eta_30",25, -2.5, 2.5,";#it{#eta^{j#gamma#gamma}");
		HStore()->createTH1F("jgamgam1_y_30",25, -2.5, 2.5,";#it{y^{j#gamma#gamma}");

		HStore()->createTH1F("jh_dR_30", 100, 0., 10.,";#it{dr}^{#gamma{1}j}");
		HStore()->createTH1F("jh_dphi_30",64, -8., 8.,";#it{dr}^{#gamma{1}j}");
		HStore()->createTH1F("jh_dy_30", 50, -5., 5.,";#it{dr}^{#gamma{1}j}");





    }
    if(_istruth){	///Initializing the hisots
		//truth photons
		HStore()->createTH1F("htrue_m",60, 100, 160,";#it{m}_{#gamma#gamma}^{truth} [GeV]");
		HStore()->createTH1F("htrue_pt",100, 0., 200.,";#it{p}^{#gamma#gamma,truth}_{T} [GeV]");
		HStore()->createTH1F("htrue_y",25, -2.5, 2.5,";#it{y}^{#gamma#gamma,truth}");
		HStore()->createTH1F("htrue_dr",100, 0., 10.,";#it{dr}^{#gamma#gamma,truth}");
		HStore()->createTH1F("htrue_dy",50, -5., 5.,";#it{dy}^{#gamma#gamma,truth}");
		HStore()->createTH1F("htrue_phi",64, -8., 8.,";#it{dphi}^{#gamma#gamma,truth}");
		//singles
		HStore()->createTH1F("p1true_pt",100, 0., 200.,";#it{p}^{#gamma,truth}_{T} [GeV]");
		HStore()->createTH1F("p2true_pt",100, 0., 200.,";#it{p}^{#gamma,truth}_{T} [GeV]");
		HStore()->createTH1F("p1true_eta",50,-2.5, 2.5,";#it{#eta}^{#gamma,truth}");
		HStore()->createTH1F("p2true_eta",50,-2.5, 2.5,";#it{#eta}^{#gamma,truth}");
		HStore()->createTH1F("p1true_phi",30,-3.2, 3.2,";#it{#phi}^{#gamma,truth}");
		HStore()->createTH1F("p2true_phi",30,-3.2, 3.2,";#it{#phi}^{#gamma,truth}");

		//fiducial
		HStore()->createTH1F("htruefid_m",60, 100, 160,";#it{m}_{#gamma#gamma}^{truth,fid} [GeV]");
		HStore()->createTH1F("htruefid_pt",100, 0., 200.,";#it{p}^{#gamma#gamma,truth,fid}_{T} [GeV]");
		HStore()->createTH1F("htruefid_y",25, -2.5, 2.5,";#it{y}^{#gamma#gamma,truth,fid}");
		HStore()->createTH1F("htruefid_dr",100, 0., 10.,";#it{dr}^{#gamma#gamma,truth,fid}");
		HStore()->createTH1F("htruefid_dy",50, -5., 5.,";#it{dy}^{#gamma#gamma,truth,fid}");
		HStore()->createTH1F("htruefid_phi",64, -8., 8.,";#it{dphi}^{#gamma#gamma,truth,fid}");
		//fiducial singles
		HStore()->createTH1F("p1truefid_pt",100, 0., 200.,";#it{p}^{#gamma,truth,fid}_{T} [GeV]");
		HStore()->createTH1F("p2truefid_pt",100, 0., 200.,";#it{p}^{#gamma,truth,fid}_{T} [GeV]");
		HStore()->createTH1F("p1truefid_eta",50,-2.5, 2.5,";#it{#eta}^{#gamma,truth,fid}");
		HStore()->createTH1F("p2truefid_eta",50,-2.5, 2.5,";#it{#eta}^{#gamma,truth,fid}");
		HStore()->createTH1F("p1truefid_phi",30,-3.2, 3.2,";#it{#phi}^{#gamma,truth,fid}");
		HStore()->createTH1F("p2truefid_phi",30,-3.2, 3.2,";#it{#phi}^{#gamma,truth,fid}");


		//fiducial with c jet
		HStore()->createTH1F("chtruefid_m",60, 100, 160,";#it{m}_{#gamma#gamma}^{truth,fid} [GeV]");
		HStore()->createTH1F("chtruefid_pt",100, 0., 200.,";#it{p}^{#gamma#gamma,truth,fid}_{T} [GeV]");
		HStore()->createTH1F("chtruefid_y",25, -2.5, 2.5,";#it{y}^{#gamma#gamma,truth,fid}");
		HStore()->createTH1F("chtruefid_dr",100, 0., 10.,";#it{dr}^{#gamma#gamma,truth,fid}");
		HStore()->createTH1F("chtruefid_dy",50, -5., 5.,";#it{dy}^{#gamma#gamma,truth,fid}");
		HStore()->createTH1F("chtruefid_phi",64, -8., 8.,";#it{dphi}^{#gamma#gamma,truth,fid}");
		//fiducial singles
		HStore()->createTH1F("cp1truefid_pt",100, 0., 200.,";#it{p}^{#gamma,truth,fid}_{T} [GeV]");
		HStore()->createTH1F("cp2truefid_pt",100, 0., 200.,";#it{p}^{#gamma,truth,fid}_{T} [GeV]");
		HStore()->createTH1F("cp1truefid_eta",50,-2.5, 2.5,";#it{#eta}^{#gamma,truth,fid}");
		HStore()->createTH1F("cp2truefid_eta",50,-2.5, 2.5,";#it{#eta}^{#gamma,truth,fid}");
		HStore()->createTH1F("cp1truefid_phi",30,-3.2, 3.2,";#it{#phi}^{#gamma,truth,fid}");
		HStore()->createTH1F("cp2truefid_phi",30,-3.2, 3.2,";#it{#phi}^{#gamma,truth,fid}");
	}




	// general - optimizations
	HStore()->createTH1F("m_yy_j25", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nojet", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nojet_inside", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//jvt weight applied
	HStore()->createTH1F("m_yy_1j_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_j20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//after selecting leading jet
// 	HStore()->createTH1F("m_yy_leadj25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_leadj20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_leadj30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_leadj35_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_dR_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("weights_after_ptetacut", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("total_weights_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");



	HStore()->createTH1F("N_j_all", 20, 0,20 ,";N_{j}");
	HStore()->createTH1F("N_j30", 20, 0,20 ,";N_{j}");

	HStore()->createTH1F("N_cjets_weighted", 20, -0.5,19 ,";N_{c}");
	HStore()->createTH1F("N_cjets", 20, -0.5,19 ,";N_{c}");

	HStore()->createTH1F("m_yy_ctag_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nonctag_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

// 	HStore()->createTH1F("NN_vs_PV",2,0,2,";matched vertices");
// 	HStore()->createTH1F("NN_vs_TRUE",2,0,2,";matched vertices");
// 	HStore()->createTH1F("PV_vs_TRUE",2,0,2,";matched vertices");
//
// 	HStore()->createTH1F("NN_vs_PV_nonctag",2,0,2,";matched vertices");
// 	HStore()->createTH1F("NN_vs_PV_ctagged",2,0,2,";matched vertices");
//
// 	HStore()->createTH1F("NN_vs_TRUE_nonctag",2,0,2,";matched vertices");
// 	HStore()->createTH1F("NN_vs_TRUE_ctagged",2,0,2,";matched vertices");
//
// 	HStore()->createTH1F("PV_vs_TRUE_nonctag",2,0,2,";matched vertices");
// 	HStore()->createTH1F("PV_vs_TRUE_ctagged",2,0,2,";matched vertices");
//
// 	HStore()->createTH1F("m_yy_pv", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_nn", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_true", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


//
// 	HStore()->createTH1F("customlead_jet_pt",400,0.,400.,";#it{p^{lead jet}_{T}");
// 	HStore()->createTH1F("pflowlead_jet_pt",400,0.,400.,";#it{p^{lead jet}_{T}");


	HStore()->createTH1F("truejets_pt",400,0.,400.,";#it{p_{T}^{j}}");
	HStore()->createTH1F("truejets_eta",50,-4,4,";#it{p_{#eta}^{j}}");
	HStore()->createTH1F("trueCjets_pt",400,0.,400.,";#it{p_{T}^{Cj}}");
	HStore()->createTH1F("trueCjets_eta",50,-4,4,";#it{p_{#eta}^{Cj}}");



/*
	HStore()->createTH1F("m_yy_dR05_mjhdefc_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjhdefc_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjhdefc_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	//with c jets
	HStore()->createTH1F("m_yy_dR05_mjh150c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh150c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165c_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR05_mjhdef_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjhdef_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjhdef_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjhdef_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjhdef_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjhdef_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjhdef_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjhdef_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjhdef_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");



    //no c jets
	HStore()->createTH1F("m_yy_dR05_mjh150_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh150_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt20_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh150_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh150_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt25_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh150_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh150_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt30_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
*/


	//HStore()->createTH1F("m_yy_tight", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_loose", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_looser", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_ofic", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	//HStore()->createTH1F("m_yy_ofic_jvt", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
		HStore()->createTH1F("m_yy_ofic_jvt_ctag", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_ofic_jvt_light", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_full_nonctag", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_bkg", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_bkg30", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_bkg35", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_bkg40", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_bkg50", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_sig30", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_sig35", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_sig40", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_sig50", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_ctag_bkg30", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_bkg35", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_bkg40", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_bkg50", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_ctag_sig30", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig35", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig40", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig50", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_sig_1j", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_sig_1jplus", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_ctag_sig_1j", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig_1jplus", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_sig", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_bkg", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_full_nonctag_jvt", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_jvt", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_bkg_jvt", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_sig_jvt", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_bkg_jvt", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig_jvt", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_full_nonctag_light", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_light", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_full_nonctag_bkg_light", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_nonctag_sig_light", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_bkg_light", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_full_ctag_sig_light", 110, 105, 160,";#it{m}_{#gamma#gamma} [GeV]");



	// HStore()->createTH1F("m_yy_full_nonctag_pv", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_full_ctag_pv", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_ofic_2taggedCtag", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_ofic_leadCtag", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_ofic_subleadCtag", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");



	//HStore()->createTH1F("m_yy_NOloose", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_NOlooser", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_NOofic", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_tight_nob", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	//	HStore()->createTH1F("m_yy_NOofic_jvt", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_NOofic_jvt_ctag", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    //HStore()->createTH1F("m_yy_NOofic_jvt_light", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

// 	HStore()->createTH1F("m_yy_tight_j30", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_loose_j30", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_looser_j30", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_ofic_j30", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
//
// 	HStore()->createTH1F("m_yy_NOloose_j30", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_NOlooser_j30", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_NOofic_j30", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
//
//
// 	HStore()->createTH1F("m_yy_tight_nob_j30", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
//
// 	HStore()->createTH1F("m_yy_allj_25gev", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_allj_30gev", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("m_yy_allj_50gev", 120, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
//
//
//
// 	HStore()->createTH1F("lead_jet_20gev",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("lead_cjet_20gev",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("lead_bjet_20gev",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("lead_ljet_20gev",150, 0, 300,"p_{T} [GeV]");
//
// 	HStore()->createTH1F("alljet_pt_25",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("alljet_pt_30",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("alljet_pt_50",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("allcjet_pt",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("allbjet_pt",150, 0, 300,"p_{T} [GeV]");
// 	HStore()->createTH1F("allljet_pt",150, 0, 300,"p_{T} [GeV]");

	//HStore()->createTH1F("lead_jet_flavour_20gev",20,-0.5,19.5,"Leadj_flavour");
	//HStore()->createTH1F("lead_jet_flavour_25gev",20,-0.5,19.5,"Leadj_flavour");
	//HStore()->createTH1F("lead_jet_flavour_30gev",20,-0.5,19.5,"Lead_flavour");


	/*
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR05_mjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR05_mjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR05_mjhdef_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR15_mjhdef_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR15_mjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR15_mjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR15_mjhdef_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR1_mjhdef_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR1_mjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

    HStore()->createTH1F("m_yy_dR1_mjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
    HStore()->createTH1F("m_yy_dR1_mjhdef_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");



    //no c jets
	HStore()->createTH1F("m_yy_dR05_mjh150_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh150_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//HStore()->createTH1F("m_yy_dR05_mjh150_pt25_ctagt60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh150_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh150_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR05_mjh155_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh155_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh155_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR05_mjh160_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh160_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh160_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh165_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh165_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR05_mjh165_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR1_mjh150_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh150_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh150_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR1_mjh155_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh155_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh155_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR1_mjh160_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh160_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh160_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh165_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh165_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_mjh165_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	//no c jets
	HStore()->createTH1F("m_yy_dR15_mjh150_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR15_mjh155_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh155_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh155_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR15_mjh160_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh160_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh160_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh165_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh165_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh165_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");



    HStore()->createTH1F("m_yy_dR05_mjh150_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh155_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh160_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_mjh165_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR1_mjh150_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh150_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh155_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh160_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_mjh165_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_mjh150_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh150_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh155_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh160_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_mjh165_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt20_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt20_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt20_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt20_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt25_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_nodR_nomjhdef_pt30_ctag60_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");



	HStore()->createTH1F("m_yy_dR05_nomjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_nomjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_nomjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_nomjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_nomjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_nomjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_nomjhdef_pt25_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_nomjhdef_pt25_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_nomjhdef_pt25_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");


	HStore()->createTH1F("m_yy_dR05_nomjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_nomjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR05_nomjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR1_nomjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_nomjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR1_nomjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

	HStore()->createTH1F("m_yy_dR15_nomjhdef_pt30_ctagtightb70_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_nomjhdef_pt30_ctagtight_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	HStore()->createTH1F("m_yy_dR15_nomjhdef_pt30_ctagloose_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

*/
//
// 	HStore()->createTH1F("Jet_index_1ctagLoose_j25", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagLooser_j25", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagOfic_j25", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagTight_j25", 25, -4.5,20.5,";j_index");
//
// 	HStore()->createTH1F("Jet_index_1ctagLoose_j30", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagLooser_j30", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagOfic_j30", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagTight_j30", 25, -4.5,20.5,";j_index");
//
// 	HStore()->createTH1F("Jet_index_1ctagLoose_j50", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_1ctagTight_j50", 25, -4.5,20.5,";j_index");
//
//
// 	HStore()->createTH1F("TruthLabel_FirstCtagLoose_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");
// 	HStore()->createTH1F("TruthLabel_FirstCtagTight_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");
//
// 	HStore()->createTH1F("TruthLabel_SecondCtagTight_j25",20,-0.5,19.5,"SecondCtagLoose_flavour");
//
 	HStore()->createTH1F("TruthLabel_NonCtagTight_j25",20,-0.5,19.5,"Nonctag_flavour");
 	HStore()->createTH1F("TruthLabel_CtagTight_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");

 	HStore()->createTH1F("N_Ctagged",20,-0.5,19.5,"Ctagged");
 	HStore()->createTH1F("N_NoNctagged",20,-0.5,19.5,"NoNctagged");
    //HStore()->createTH1F("TruthLabel_NonCtagTight_lead",20,-0.5,19.5,"Nonctag_flavour");
    //HStore()->createTH1F("TruthLabel_NonCtagTight_sublead",20,-0.5,19.5,"Nonctag_flavour");
//
// 	HStore()->createTH1F("TruthLabel_FirstCtagLoose_j30",20,-0.5,19.5,"FirstCtagLoose_flavour");
// 	HStore()->createTH1F("TruthLabel_FirstCtagTight_j30",20,-0.5,19.5,"FirstCtagLoose_flavour");
//
// 	HStore()->createTH1F("TruthLabel_LeadCtagLoose_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");
// 	HStore()->createTH1F("TruthLabel_LeadCtagTight_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");
//
// 	HStore()->createTH1F("TruthLabel_LeadCtagLoose_j30",20,-0.5,19.5,"FirstCtagLoose_flavour");
// 	HStore()->createTH1F("TruthLabel_LeadCtagTight_j30",20,-0.5,19.5,"FirstCtagLoose_flavour");
//
// 	HStore()->createTH1F("TruthLabel_AllCtagLoose_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");
// 	HStore()->createTH1F("TruthLabel_AllCtagTight_j25",20,-0.5,19.5,"FirstCtagLoose_flavour");
//
// 	HStore()->createTH1F("TruthLabel_AllCtagLoose_j30",20,-0.5,19.5,"FirstCtagLoose_flavour");
// 	HStore()->createTH1F("TruthLabel_AllCtagTight_j30",20,-0.5,19.5,"FirstCtagLoose_flavour");
//
//
// 	HStore()->createTH1F("dR_y1_j", 100, 0., 10.,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("dR_y2_j", 100, 0., 10.,";#it{m}_{#gamma#gamma} [GeV]");
//
// 	HStore()->createTH1F("p1pt", 100, 0., 200.,";#it{m}_{#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("p2pt", 100, 0., 200.,";#it{m}_{#gamma#gamma} [GeV]");

    //HStore()->createTH1F("lead_jet_pt", 100, 0, 200,";Lead jet p_{T}[GeV]");
//
//     HStore()->createTH1F("m_yy_HighPt_jvt",60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
//     HStore()->createTH1F("m_yy_JgamgamM_jvt",60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
//
// 	HStore()->createTH1F("leadCjet_pt",300, 0, 300,";LeadC jet p_{T}[GeV]");
// 	HStore()->createTH1F("leadCjet_eta",25, -2.5, 2.5,";#it{#eta}^{j}");
// 	HStore()->createTH1F("leadCjet_y",25, -2.5, 2.5,";#it{y}^{j}");
// 	HStore()->createTH1F("leadCjet_phi",30,-3.2, 3.2,";#it{#phi}^{j}");
// 	HStore()->createTH1F("leadCjet_m",400, 0, 400,";#it{m}_{j} [GeV]");
//
// 	HStore()->createTH1F("jgamgam_m",200, 0, 400,";#it{m}_{j#gamma#gamma} [GeV]");
// 	HStore()->createTH1F("jgamgam_pt",200, 0, 400,";#it{p}^{j#gamma#gamma}_{T} [GeV]");
// 	HStore()->createTH1F("jgamgam_eta",25, -2.5, 2.5,";#it{#eta}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_y",25, -2.5, 2.5,";#it{y}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_phi",30,-3.2, 3.2,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_dR",100, 0., 10.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_dphi",64, -8., 8.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_dy",50, -5., 5.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_dR_b",100, 0., 10.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_dphi_b",64, -8., 8.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_dy_b",50, -5., 5.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("jgamgam_m_b",400, 0, 400,";#it{#phi}^{j#gamma#gamma}");


	//after passing no b tagged jets
//	HStore()->createTH1F("m_yy_j25_mv260_jvt", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");
	//after passing the c tagging - noefficiency applied
//	HStore()->createTH1F("m_yy_j25_mv260_ctag", 60, 100, 160,";#it{m}_{#gamma#gamma} [GeV]");

//
//
// 	HStore()->createTH1F("dR_j1y1",100, 0., 10.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("dR_j1y2",100, 0., 10.,";#it{#phi}^{j#gamma#gamma}");
//
//
// 	HStore()->createTH1F("dR_j1y1_30",100, 0., 10.,";#it{#phi}^{j#gamma#gamma}");
// 	HStore()->createTH1F("dR_j1y2_30",100, 0., 10.,";#it{#phi}^{j#gamma#gamma}");
//
//
// 	HStore()->createTH1F("Jet_index_drcounts", 25, -4.5,20.5,";j_index");
// 	HStore()->createTH1F("Jet_index_drcounts_30", 25, -4.5,20.5,";j_index");



//
//
//
//     //the same but with applying the efficiency as event weight
// 	//other - weights for the events
	HStore()->createTH1F("total_weights_calc", 100,-1,2,";weights_calc_tot");
    HStore()->createTH1F("total_weights_Ntup",100,-1,2,";weights_ntup_tot");
    HStore()->createTH1F("jvt_plus_total",100,-1,2,";weights_plus_jvt");
    HStore()->createTH1F("jvt_weight",1000,-1,2,";jvt_weight");

//
//     HStore()->createTH1F("truth_tagged_C",2,-0.5,1.5,";Lead_C");
//     HStore()->createTH1F("truth_tagged_B",2,-0.5,1.5,";Lead_B");
//     HStore()->createTH1F("truth_tagged_L",2,-0.5,1.5,";Lead_L");
//
//     HStore()->createTH1F("truth_jets_AfterCtag_c",4,0,4,";Truth Label");
//     HStore()->createTH1F("truth_jets_AfterCtag_b",4,0,4,";Truth Label");
//     HStore()->createTH1F("truth_jets_AfterCtag_l",4,0,4,";Truth Label");
//

	//From the hist store
	for (auto *histo : m_hStore->getListOfHistograms())
		{wk()->addOutput(histo);}




    if(_Dooptimization){

		for(int k = 0; k<bins; k++){
			for(int j =0; j<fbins; j++){
				masses_mv60[k][j] = new TH1D(Form("masses_mv60_%d_%d",k,j),"masses_mv60",60, 100, 160);
				masses_mv60nonc[k][j] = new TH1D(Form("masses_mv60nonc_%d_%d",k,j),"masses_mv60nonc",60, 100, 160);
				// masses_mv70[k][j] = new TH1D(Form("masses_mv70_%d_%d",k,j),"masses_mv70",60, 100, 160);
// 				masses_mv77[k][j] = new TH1D(Form("masses_mv77_%d_%d",k,j),"masses_mv77",60, 100, 160);
//
// 				//for the efficiency
// 				masses_mv60_nodl1[k][j] = new TH1D(Form("masses_mv60_nodl1_%d_%d",k,j),"masses_mv60_nodl1",60, 100, 160);
//
//
// 				//jvt weight
// 				masses_mv60_jvt[k][j] = new TH1D(Form("masses_mv60_jvt_%d_%d",k,j),"masses_mvjvt60",60, 100, 160);
// 				masses_mv70_jvt[k][j] = new TH1D(Form("masses_mv70_jvt_%d_%d",k,j),"masses_mvjvt70",60, 100, 160);
// 				masses_mv77_jvt[k][j] = new TH1D(Form("masses_mv77_jvt_%d_%d",k,j),"masses_mvjvt77",60, 100, 160);
// 				//weighted with fractions
// 				masses_mv60w[k][j] = new TH1D(Form("masses_mv60w_%d_%d",k,j),"masses_mv60w",60, 100, 160);
// 				masses_mv70w[k][j] = new TH1D(Form("masses_mv70w_%d_%d",k,j),"masses_mv70w",60, 100, 160);
// 				masses_mv77w[k][j] = new TH1D(Form("masses_mv77w_%d_%d",k,j),"masses_mv77w",60, 100, 160);
//
//  				//efficiencies and numbers of jets
// 				N_recob[k][j] = new TH1D(Form("N_recob_%d_%d",k,j),"N_recob",20, -0.5, 19.5);
// 				N_recoc[k][j] = new TH1D(Form("N_recoc_%d_%d",k,j),"N_recoc",20, -0.5, 19.5);
// 				N_recol[k][j] = new TH1D(Form("N_recol_%d_%d",k,j),"N_recol",20, -0.5, 19.5);
// 				/////
// 				N_tagb[k][j]  = new TH1D(Form("N_tagb_%d_%d",k,j),"N_tagb",20, -0.5, 19.5);
// 				N_tagc[k][j]  = new TH1D(Form("N_tagc_%d_%d",k,j),"N_tagc",20, -0.5, 19.5);
// 				N_tagl[k][j]  = new TH1D(Form("N_tagl_%d_%d",k,j),"N_tagl",20, -0.5, 19.5);
// 				//////
// 				eff_bmiss[k][j]  = new TH1D(Form("eff_bmiss_%d_%d",k,j),"eff_bmiss",20, -0.05, 1.05);
// 				eff_lmiss[k][j]	 = new TH1D(Form("eff_lmiss_%d_%d",k,j),"eff_lmiss",20, -0.05, 1.05);
// 				eff_cmiss[k][j]  = new TH1D(Form("eff_cmiss_%d_%d",k,j),"eff_cmiss",20, -0.05, 1.05);
// 				//
// 				eff_bmissw[k][j]  = new TH1D(Form("eff_bmissw_%d_%d",k,j),"eff_bmissw",20, -0.05, 1.05);
// 				eff_lmissw[k][j]  = new TH1D(Form("eff_lmissw_%d_%d",k,j),"eff_lmissw",20, -0.05, 1.05);
// 				eff_cmissw[k][j]  = new TH1D(Form("eff_cmissw_%d_%d",k,j),"eff_cmissw",20, -0.05, 1.05);
// 				//
// 				eff_bmiss_nosf[k][j]  = new TH1D(Form("eff_bmiss_nosf_%d_%d",k,j),"eff_bmiss_nosf",20, -0.05, 1.05);
// 				eff_lmiss_nosf[k][j]  = new TH1D(Form("eff_lmiss_nosf_%d_%d",k,j),"eff_lmiss_nosf",20, -0.05, 1.05);
// 				eff_cmiss_nosf[k][j]  = new TH1D(Form("eff_cmiss_nosf_%d_%d",k,j),"eff_cmiss_nosf",20, -0.05, 1.05);
// 				//
// 				eff_bmiss_init[k][j]  = new TH1D(Form("eff_bmiss_init_%d_%d",k,j),"eff_bmiss_init",20, -0.05, 1.05);
// 				eff_lmiss_init[k][j]  = new TH1D(Form("eff_lmiss_init_%d_%d",k,j),"eff_lmiss_init",20, -0.05, 1.05);
// 				eff_cmiss_init[k][j]  = new TH1D(Form("eff_cmiss_init_%d_%d",k,j),"eff_cmiss_init",20, -0.05, 1.05);



				//////-----------------------------------------
				wk()->addOutput(masses_mv60[k][j]);
				wk()->addOutput(masses_mv60nonc[k][j]);
				// wk()->addOutput(masses_mv70[k][j]);
// 				wk()->addOutput(masses_mv77[k][j]);
//
// 				wk()->addOutput(masses_mv60_nodl1[k][j]);
// 				//jvt
// 				wk()->addOutput(masses_mv60_jvt[k][j]);
// 				wk()->addOutput(masses_mv70_jvt[k][j]);
// 				wk()->addOutput(masses_mv77_jvt[k][j]);
// 				//weights frac
// 				wk()->addOutput(masses_mv60w[k][j]);
// 				wk()->addOutput(masses_mv70w[k][j]);
// 				wk()->addOutput(masses_mv77w[k][j]);
// 				//fracs
// 				wk()->addOutput(N_recob[k][j]);
// 				wk()->addOutput(N_recoc[k][j]);
// 				wk()->addOutput(N_recol[k][j]);
// 				///
// 				wk()->addOutput(N_tagb[k][j]);
// 				wk()->addOutput(N_tagc[k][j]);
// 				wk()->addOutput(N_tagl[k][j]);
// 				///
// 				wk()->addOutput(eff_bmiss[k][j]);
// 				wk()->addOutput(eff_cmiss[k][j]);
// 				wk()->addOutput(eff_lmiss[k][j]);
// 				///
// 				wk()->addOutput(eff_bmissw[k][j]);
// 				wk()->addOutput(eff_cmissw[k][j]);
// 				wk()->addOutput(eff_lmissw[k][j]);
// 				//
// 				wk()->addOutput(eff_bmiss_init[k][j]);
// 				wk()->addOutput(eff_lmiss_init[k][j]);
// 				wk()->addOutput(eff_cmiss_init[k][j]);
// 				//
// 				wk()->addOutput(eff_bmiss_nosf[k][j]);
// 				wk()->addOutput(eff_lmiss_nosf[k][j]);
// 				wk()->addOutput(eff_cmiss_nosf[k][j]);




				//if(k==0 and (j%4==0)){//for 41 bin
				//if(k==0 and (j%10==0)){
				//   	dl1_60[j] = new TH1D(Form("dl1_60_%d",j),"dl1_60",101,-5.5,5.5);
					// dl1_70[j] = new TH1D(Form("dl1_70_%d",j),"dl1_70",101,-5.5,5.5);
// 					dl1_77[j] = new TH1D(Form("dl1_77_%d",j),"dl1_77",101,-5.5,5.5);
//
// 					//fracs
// 					dl1_60w[j] = new TH1D(Form("dl1_60w_%d",j),"dl1_60w",101,-5.5,5.5);
// 					dl1_70w[j] = new TH1D(Form("dl1_70w_%d",j),"dl1_70w",101,-5.5,5.5);
// 					dl1_77w[j] = new TH1D(Form("dl1_77w_%d",j),"dl1_77w",101,-5.5,5.5);

					//wk()->addOutput(dl1_60[j]);
					// wk()->addOutput(dl1_70[j]);
// 					wk()->addOutput(dl1_77[j]);
//
// 					//fracs
// 					wk()->addOutput(dl1_60w[j]);
// 					wk()->addOutput(dl1_70w[j]);
// 					wk()->addOutput(dl1_77w[j]);

				//}

			}

		}

/*
		//efficiencies 2d
		j_ctrue = new TH2D("j_ctrue","j_ctrue",100,-5.,5.,100,0.,1.);
		j_btrue = new TH2D("j_btrue","j_btrue",100,-5.,5.,100,0.,1.);
		j_ltrue = new TH2D("j_ltrue","j_ltrue",100,-5.,5.,100,0.,1.);

		c_ctrue = new TH2D("c_ctrue","c_ctrue",100,-5.,5.,100,0.,1.);
		c_btrue = new TH2D("c_btrue","c_btrue",100,-5.,5.,100,0.,1.);
		c_ltrue = new TH2D("c_ltrue","c_ltrue",100,-5.,5.,100,0.,1.);

		wk()->addOutput(j_ctrue);
		wk()->addOutput(j_btrue);
		wk()->addOutput(j_ltrue);

		wk()->addOutput(c_ctrue);
		wk()->addOutput(c_btrue);
		wk()->addOutput(c_ltrue);
*/
 	}//do the optimization



	return EL::StatusCode::SUCCESS;
}

//checking reco plots
EL::StatusCode HcAnalyse :: doreco(const xAOD::PhotonContainer *photons,const xAOD::JetContainer *jets)
{
	auto gamlead = (*photons)[0];
	auto gamsublead = (*photons)[1];
	//1 photon //2 photon
	float y1_pt  = gamlead->pt();   float y2_pt  = gamsublead->pt();
	float y1_eta = gamlead->eta();  float y2_eta = gamsublead->eta();
	float y1_phi = gamlead->phi();  float y2_phi = gamsublead->phi();

	//Higgs candidates
	TLorentzVector ph1 = gamlead->p4(), ph2 = gamsublead->p4();
	ph1 *= HC::invGeV; ph2 *= HC::invGeV; TLorentzVector hcand = ph1+ph2;


	auto leadjet = (*jets)[0];
	TLorentzVector jet = leadjet->p4();
	jet *= HC::invGeV;

	TLorentzVector jgamgam1 = jet + hcand;

	//if(jet.Pt()>25 && fabs(jet.PseudoRapidity())<2.5){
	if(jet.Pt()>25){
	//higgs candates
	HStore()->fillTH1F("myy",hcand.M(),weight_analyse);
	HStore()->fillTH1F("ptyy",hcand.Pt(),weight_analyse);
	HStore()->fillTH1F("yyy",hcand.Rapidity(),weight_analyse);
	HStore()->fillTH1F("dryy",ph1.DeltaR(ph2),weight_analyse);
	HStore()->fillTH1F("dyyy",ph1.Rapidity() - ph2.Rapidity(),weight_analyse);
	HStore()->fillTH1F("dphiyy",ph1.DeltaPhi(ph2),weight_analyse);
	//photons
	HStore()->fillTH1F("y1pt",y1_pt * HC::invGeV,weight_analyse);
	HStore()->fillTH1F("y2pt",y2_pt * HC::invGeV,weight_analyse);
	HStore()->fillTH1F("y1eta",y1_eta,weight_analyse);
	HStore()->fillTH1F("y2eta",y2_eta,weight_analyse);
	HStore()->fillTH1F("y1phi",y1_phi,weight_analyse);
	HStore()->fillTH1F("y2phi",y2_phi,weight_analyse);
	//Drs
	HStore()->fillTH1F("dRjgam1",jet.DeltaR(ph1),weight_analyse);
	HStore()->fillTH1F("dRjgam2",jet.DeltaR(ph2),weight_analyse);
	HStore()->fillTH1F("dphijgam1",jet.DeltaPhi(ph1),weight_analyse);
	HStore()->fillTH1F("dphijgam2",jet.DeltaPhi(ph2),weight_analyse);
	//jet kinematics
	HStore()->fillTH1F("jpt",jet.Pt(),weight_analyse);
	HStore()->fillTH1F("jy",jet.Rapidity(),weight_analyse);
	HStore()->fillTH1F("jeta",jet.PseudoRapidity(),weight_analyse);
	//j+gam+gam
	HStore()->fillTH1F("jgamgam1_pt",jgamgam1.Pt(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_m",jgamgam1.M(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_phi",jgamgam1.Phi(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_eta",jgamgam1.PseudoRapidity(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_y",jgamgam1.Rapidity(),weight_analyse);
	HStore()->fillTH1F("jh_dR",jet.DeltaR(hcand),weight_analyse);
	HStore()->fillTH1F("jh_dphi",jet.DeltaPhi(hcand),weight_analyse);
	HStore()->fillTH1F("jh_dy",jet.Rapidity()-hcand.Rapidity(),weight_analyse);
    }


	if(jet.Pt()>30 && fabs(jet.PseudoRapidity())<2.5){

	TLorentzVector jgamgam1 = jet + hcand;

	HStore()->fillTH1F("dRjgam1_30",jet.DeltaR(ph1),weight_analyse);
	HStore()->fillTH1F("dRjgam2_30",jet.DeltaR(ph2),weight_analyse);
	HStore()->fillTH1F("dphijgam1_30",jet.DeltaPhi(ph1),weight_analyse);
	HStore()->fillTH1F("dphijgam2_30",jet.DeltaPhi(ph2),weight_analyse);
	//jet kinematics
	HStore()->fillTH1F("jeta_30",jet.PseudoRapidity(),weight_analyse);
	//j+gam+gam
	HStore()->fillTH1F("jgamgam1_pt_30",jgamgam1.Pt(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_m_30",jgamgam1.M(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_phi_30",jgamgam1.Phi(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_eta_30",jgamgam1.PseudoRapidity(),weight_analyse);
	HStore()->fillTH1F("jgamgam1_y_30",jgamgam1.Rapidity(),weight_analyse);
	HStore()->fillTH1F("jh_dR_30",jet.DeltaR(hcand),weight_analyse);
	HStore()->fillTH1F("jh_dphi_30",jet.DeltaPhi(hcand),weight_analyse);
	HStore()->fillTH1F("jh_dy_30",jet.Rapidity()-hcand.Rapidity(),weight_analyse);

	}





	return EL::StatusCode::SUCCESS;

}

//--------------------------------------------------------------------------------------------------
//checking truth plots
EL::StatusCode HcAnalyse :: dotruth(const xAOD::TruthParticleContainer *truphotons,const xAOD::EventInfo* HgamTrutheventInfo)
{

	//checking the number of photons
	if (truphotons->size() > 1){
		auto lead = (*truphotons)[0];
		auto sublead = (*truphotons)[1];
		//1 photon //2 photon
		float ph1_pt  = lead->pt();		float ph2_pt  = sublead->pt();
		float ph1_eta = lead->eta();	float ph2_eta = sublead->eta();
		float ph1_phi = lead->phi();	float ph2_phi = sublead->phi();

		//Higgs candidates
		TLorentzVector p1 = lead->p4(), p2 = sublead->p4();
		p1 *= HC::invGeV; p2 *= HC::invGeV; TLorentzVector h = p1+p2;

		//higgs truth
		HStore()->fillTH1F("htrue_m",h.M());
		HStore()->fillTH1F("htrue_pt",h.Pt());
		HStore()->fillTH1F("htrue_y",h.Rapidity());
		HStore()->fillTH1F("htrue_dr",p1.DeltaR(p2));
		HStore()->fillTH1F("htrue_dy",p1.Rapidity() - p2.Rapidity());
		HStore()->fillTH1F("htrue_phi",p1.DeltaPhi(p2));
		//photons truth
		HStore()->fillTH1F("p1true_pt",ph1_pt * HC::invGeV);
		HStore()->fillTH1F("p2true_pt",ph2_pt * HC::invGeV);
		HStore()->fillTH1F("p1true_eta",ph1_eta);
		HStore()->fillTH1F("p2true_eta",ph2_eta);
		HStore()->fillTH1F("p1true_phi",ph1_phi);
		HStore()->fillTH1F("p2true_phi",ph2_phi);

		if(HgamTrutheventInfo->auxdataConst<char>("isFiducial")>0){
			//1 photons
			float ph1pt  = lead->pt();	float ph2pt  = sublead->pt();
			float ph1eta = lead->eta();	float ph2eta = sublead->eta();
			float ph1phi = lead->phi();	float ph2phi = sublead->phi();

			//Higgs candidates
			TLorentzVector p1fid = lead->p4(), p2fid = sublead->p4();
			p1fid *= HC::invGeV; p2fid *= HC::invGeV; TLorentzVector hfid = p1fid+p2fid;

			//fiducial truth
			HStore()->fillTH1F("htruefid_m",hfid.M());
			HStore()->fillTH1F("htruefid_pt",hfid.Pt());
			HStore()->fillTH1F("htruefid_y",hfid.Rapidity());
			HStore()->fillTH1F("htruefid_dr",p1fid.DeltaR(p2fid));
			HStore()->fillTH1F("htruefid_dy",p1fid.Rapidity() - p2fid.Rapidity());
			HStore()->fillTH1F("htruefid_phi",p1fid.DeltaPhi(p2fid));
			HStore()->fillTH1F("p1truefid_pt",ph1pt * HC::invGeV);
			HStore()->fillTH1F("p2truefid_pt",ph2pt * HC::invGeV);
			HStore()->fillTH1F("p1truefid_eta",ph1eta);
			HStore()->fillTH1F("p2truefid_eta",ph2eta);
			HStore()->fillTH1F("p1truefid_phi",ph1phi);
			HStore()->fillTH1F("p2truefid_phi",ph2phi);

			if(HgamTrutheventInfo->auxdataConst<int>("N_j_ctag25")==1){
			  	//1 photons
				float cph1pt  = lead->pt();		float cph2pt  = sublead->pt();
				float cph1eta = lead->eta();	float cph2eta = sublead->eta();
				float cph1phi = lead->phi();	float cph2phi = sublead->phi();

				//Higgs candidates
				TLorentzVector cp1fid = lead->p4(), cp2fid = sublead->p4();
				cp1fid *= HC::invGeV; cp2fid *= HC::invGeV; TLorentzVector chfid = cp1fid+cp2fid;

				//fiducial truth
				HStore()->fillTH1F("chtruefid_m",chfid.M());
				HStore()->fillTH1F("chtruefid_pt",chfid.Pt());
				HStore()->fillTH1F("chtruefid_y",chfid.Rapidity());
				HStore()->fillTH1F("chtruefid_dr",cp1fid.DeltaR(cp2fid));
				HStore()->fillTH1F("chtruefid_dy",cp1fid.Rapidity() - cp2fid.Rapidity());
				HStore()->fillTH1F("chtruefid_phi",cp1fid.DeltaPhi(cp2fid));
				HStore()->fillTH1F("cp1truefid_pt",cph1pt * HC::invGeV);
				HStore()->fillTH1F("cp2truefid_pt",cph2pt * HC::invGeV);
				HStore()->fillTH1F("cp1truefid_eta",cph1eta);
				HStore()->fillTH1F("cp2truefid_eta",cph2eta);
				HStore()->fillTH1F("cp1truefid_phi",cph1phi);
				HStore()->fillTH1F("cp2truefid_phi",cph2phi);

			}//end of c tagged jet

		}//end of fiducial statement

	}//end of "if" statement

return EL::StatusCode::SUCCESS;
}

//double HcAnalyse::DeltaR(double phi1, double eta1, double phi2, double eta2)
//{
//  double deta=eta1-eta2;
//  double dphi=phi1-phi2;
//  while (dphi < -TMath::Pi() ) dphi += 2*TMath::Pi();
//  while (dphi >  TMath::Pi() ) dphi -= 2*TMath::Pi();
//  return std::sqrt(deta*deta+dphi*dphi);
//}


//if you want to use not the leading jet
//unsigned int thisjet=0;
//				for(;thisjet<jets->size();++thisjet)
//				{
//					if( ((*jets)[thisjet]->pt()*GeV > 25)  && ( fabs( (*jets)[thisjet]->eta() )<2.5 ) ) break;
//				}
//				auto leadjet = (*jets)[thisjet];
//



//Here lies the code for the truth tagging which is now irrelevant - RIP
/*
	//define the variables -----------------------------------------//
		int countcjets_b60 = 0;
		float epsilon_b = 0.0, epsilon_c = 0.0, epsilon_l = 0.0;

        if(isMC){
			// ------ doing the c tagging for the jets using the selected points dl1=0.39 and frac 0.22 ---- ////
			std::cout << "-------------------------------------------------------" << std::endl;
			std::cout << " Doing the c tagging using dl = 0.4 and frac = 0.23    " << std::endl;
			std::cout << "-------------------------------------------------------" << std::endl;

       		//define the variables -----------------------------------------//
			int countTrueBmatch=0, countTrueCmatch=0, countTrueLmatch=0;
			int bmistag =0, ctag=0, lmistag =0;
			//--------------------------------------------------------------//

			//-----------------First lets count number of true c,b, and light jets and count the Number of misidentified-----------------//
    		for(auto recoj: *jets){
    			if(fabs(recoj->eta()) > 2.5)continue;//this is after version mistag_v3,v4
    			if(recoj->pt() < 25 * HC::GeV)continue;//this is after version mistag_v3,v4
    			bool isCjet = false, isBjet = false, isLjet = false;
    			//Checking the truth labeling
    			if(recoj->auxdataConst<int>("HadronConeExclTruthLabelID") == 5){isBjet = true;}
    			if(recoj->auxdataConst<int>("HadronConeExclTruthLabelID") == 4){isCjet = true;}
    			if(recoj->auxdataConst<int>("HadronConeExclTruthLabelID") == 0){isLjet = true;}
    			//setting the booleans
		    	bool truth_b = false; bool truth_l = false; bool truth_c = false;

				if(isBjet){truth_b = true; countTrueBmatch++;}
				if(isCjet){truth_c = true; countTrueCmatch++;}
				if(isLjet){truth_l = true; countTrueLmatch++;}

		    	//DL1
				double pu_j = recoj->auxdataConst<double>("DL1_pu");
				double pc_j = recoj->auxdataConst<double>("DL1_pc");
				double pb_j = recoj->auxdataConst<double>("DL1_pb");
				//Dl
				double dl1discriminant = discriminant(pc_j,pb_j,pu_j,0.23);
				if(truth_b && (dl1discriminant > 0.4)){bmistag++;}
				if(truth_c && (dl1discriminant > 0.4)){ctag++; }
				if(truth_l && (dl1discriminant > 0.4)){lmistag++; }

			}//end of reco jet loop


			int nj = jets->size();
			HStore()->fillTH1F("Nj_total",nj);
			HStore()->fillTH1F("N_trub",countTrueBmatch);
			HStore()->fillTH1F("N_truc",countTrueCmatch);
			HStore()->fillTH1F("N_trul",countTrueLmatch);

			HStore()->fillTH1F("N_trub_missc",bmistag);
			HStore()->fillTH1F("N_truc_missc",ctag);
			HStore()->fillTH1F("N_trul_missc",lmistag);


			//-----------------Second lets make the efficiency for the misidentification-----------------//
        	//misidentification B
			if(countTrueBmatch==0) {epsilon_b = 0.0;}
			else { epsilon_b = (bmistag*1.0)/(1.0*countTrueBmatch);}

			//identification C
			if(countTrueCmatch ==0) {epsilon_c = 0.0;}
        	else { epsilon_c = (ctag*1.0)/(1.0*countTrueCmatch); }

	    	//misidentification L
			if(countTrueLmatch==0) {epsilon_l = 0.0;}
			else { epsilon_l = (lmistag*1.0)/(1.0*countTrueLmatch); }

			//filling the efficiency histograms
			if(countTrueBmatch > 0){
			HStore()->fillTH1F("eff_bmiss",epsilon_b);                      HStore()->fillTH1F("eff_bmissw",epsilon_b,weightjvt);
			HStore()->fillTH1F("eff_bmiss_nosf",epsilon_b,tot_weight_nosf); HStore()->fillTH1F("eff_bmiss_init",epsilon_b,initial_weight);}

			if(countTrueCmatch > 0){
			HStore()->fillTH1F("eff_c",epsilon_c);      					HStore()->fillTH1F("eff_cw",epsilon_c,weightjvt);
			HStore()->fillTH1F("eff_c_nosf",epsilon_c,tot_weight_nosf);  	HStore()->fillTH1F("eff_c_init",epsilon_c,initial_weight);}


			if(countTrueLmatch > 0){
			HStore()->fillTH1F("eff_lmiss",epsilon_l);   					HStore()->fillTH1F("eff_lmissw",epsilon_l,weightjvt);
			HStore()->fillTH1F("eff_lmiss_nosf",epsilon_l,tot_weight_nosf); HStore()->fillTH1F("eff_lmiss_init",epsilon_l,initial_weight);}
        }//end of isMC


           // bool isCjet_tagged = false, isBjet_tagged = false, isLjet_tagged = false;
		    //if(leadjet->auxdataConst<int>("HadronConeExclTruthLabelID") == 4) {isCjet_tagged = true;}
		   // if(leadjet->auxdataConst<int>("HadronConeExclTruthLabelID") == 5) {isBjet_tagged = true;}
		   // if(leadjet->auxdataConst<int>("HadronConeExclTruthLabelID") == 0) {isLjet_tagged = true;}
		   // HStore()->fillTH1F("truth_jets_AfterCtag_c",isCjet_tagged);
		   // HStore()->fillTH1F("truth_jets_AfterCtag_b",isBjet_tagged);
		   // HStore()->fillTH1F("truth_jets_AfterCtag_l",isLjet_tagged);



		//  bool isCjet = false, isBjet = false, isLjet = false;;
		//  if(leadjet->auxdataConst<int>("HadronConeExclTruthLabelID") == 4){ isCjet = true;} HStore()->fillTH1F("truth_tagged_C",isCjet);
        //  if(leadjet->auxdataConst<int>("HadronConeExclTruthLabelID") == 5){ isBjet = true;} HStore()->fillTH1F("truth_tagged_B",isBjet);
    	//  if(leadjet->auxdataConst<int>("HadronConeExclTruthLabelID") == 0){ isLjet = true;} HStore()->fillTH1F("truth_tagged_L",isLjet);


*/



//Something suspicious
/*
        if(mv2_60==false && dl1_discriminant > 0.4 && dr05_nomhj==true){HStore()->fillTH1F("m_yy_dR05_mjhdefc_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr10_nomhj==true){HStore()->fillTH1F("m_yy_dR1_mjhdefc_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr15_nomhj==true){HStore()->fillTH1F("m_yy_dR15_mjhdefc_jvt",mass_gev,weightjvt20);}

        if(mv2_60==false && dl1_discriminant > 0.4 && dr05_mhj150==true){HStore()->fillTH1F("m_yy_dR05_mjh150c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr05_mhj155==true){HStore()->fillTH1F("m_yy_dR05_mjh155c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr05_mhj160==true){HStore()->fillTH1F("m_yy_dR05_mjh160c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr05_mhj165==true){HStore()->fillTH1F("m_yy_dR05_mjh165c_jvt",mass_gev,weightjvt20);}

        if(mv2_60==false && dl1_discriminant > 0.4 && dr10_mhj150==true){HStore()->fillTH1F("m_yy_dR1_mjh150c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr10_mhj155==true){HStore()->fillTH1F("m_yy_dR1_mjh155c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr10_mhj160==true){HStore()->fillTH1F("m_yy_dR1_mjh160c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr10_mhj165==true){HStore()->fillTH1F("m_yy_dR1_mjh165c_jvt",mass_gev,weightjvt20);}

        if(mv2_60==false && dl1_discriminant > 0.4 && dr15_mhj150==true){HStore()->fillTH1F("m_yy_dR15_mjh150c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr15_mhj155==true){HStore()->fillTH1F("m_yy_dR15_mjh155c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr15_mhj160==true){HStore()->fillTH1F("m_yy_dR15_mjh160c_jvt",mass_gev,weightjvt20);}
        if(mv2_60==false && dl1_discriminant > 0.4 && dr15_mhj165==true){HStore()->fillTH1F("m_yy_dR15_mjh165c_jvt",mass_gev,weightjvt20);}
*/



//other code


//
//
// 		int n_btagged_jets = 0; // initialize the counter for b-tagged jets
// 		int n_nonbtagged_jets = 0; // initialize the counter for non-btagged jets
// 		int jetcount = 0;
// 		int btagged_jet_flavour = 0;
// 		int nonbtagged_jet_flavour = 0;
//
// 		for(auto recoj: *jets){
//     		if(fabs(recoj->eta()) > 2.5) continue;
//     		if(recoj->pt() < 25 * HC::GeV) continue;
//     		jetcount++;
//
//     		if(jetcount>=3)break;
//
//     		if(HC::DR(recoj,pho_lead) < _dR_cut || HC::DR(recoj,pho_sublead) < _dR_cut) continue;
//
//
//
//     		double pu_jr = recoj->auxdataConst<double>("DL1r_pu");
//     		double pc_jr = recoj->auxdataConst<double>("DL1r_pc");
//     		double pb_jr = recoj->auxdataConst<double>("DL1r_pb");
//
//     		double dl1_discriminant_DL1Ofic = discriminant(pc_jr, pb_jr, pu_jr, 0.2);
//
//     		if(dl1_discriminant_DL1Ofic > 1) {
//         		n_btagged_jets++;
//         		btagged_jet_flavour = recoj->auxdataConst<int>("HadronConeExclTruthLabelID");
//         		break;
//         	}
//     		else if(dl1_discriminant_DL1Ofic < 1 ) {
//         		n_nonbtagged_jets++;
//         		nonbtagged_jet_flavour = recoj->auxdataConst<int>("HadronConeExclTruthLabelID");
//         		break;
//     		}
//
//         }
//
//     	if(n_btagged_jets > 0 ) { HStore()->fillTH1F("m_yy_ofic_jvt_ctag",mass_gev,weightjvt_ctag); HStore()->fillTH1F("TruthLabel_CtagTight_j25", btagged_jet_flavour);}
//     	else if(n_nonbtagged_jets>0) {HStore()->fillTH1F("m_yy_NOofic_jvt_ctag",mass_gev,weightjvt_ctag); HStore()->fillTH1F("TruthLabel_NonCtagTight_j25", nonbtagged_jet_flavour); }
//
//     	HStore()->fillTH1F("N_NoNctagged",n_nonbtagged_jets);
//     	HStore()->fillTH1F("N_Ctagged",n_btagged_jets);
