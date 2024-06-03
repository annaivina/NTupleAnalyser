#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/CondorDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopGrid/GridDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include "EventLoop/OutputStream.h"
#include "EventLoopAlgs/NTupleSvc.h"

#include <TSystem.h>

#include "SampleHandler/ScanDir.h"

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>


#include "NTupAnalyser/HcAnalyse.h"


int main(int argc, char *argv[])
{

	//defining the variables
	std::string inputDir;
	std::string submitDir;
	std::string InDS;
	std::string OutDS;
	std::string output_file_name;
	std::string rootFile;
	std::string mcname;
	int Nevents;
	bool isGridJob;
	bool ispTcut;
	bool signal;
	bool doCutFlow;
	bool optimization;
	bool istruth;
	bool isreco;
	int nFilesPerJob;
	float phPtcut;
	float jetpT;
	float jeteta;
	float frac;
	float cut;
	float dR_cut;
	string reco_jet_collection;
	string b_tagger;



	//Boost configuration
	//1) command line only: options can only be given on command line, not in config file
	boost::program_options::options_description cmd_only_options("command line only options");
	std::string config_file_name;

	cmd_only_options.add_options() //note unusual syntax when adding options!
		("help,h","produce help message")
		("config,c",boost::program_options::value<std::string>(&config_file_name),"name of configuration file");

	//2) main options: most likely to be set by user, can be specified both via command line or config
	//explaination is included in help message
	boost::program_options::options_description main_options("main options");

	main_options.add_options()
		("output_file",boost::program_options::value<std::string>(&output_file_name)->default_value("ntuple"),"name of output root file")
		("inputDir",boost::program_options::value<std::string>(&inputDir)->default_value("/afs/cern.ch/user/a/aivina/cernbox/MC_testing_mass/"),"name of input directory containing all files")
		("submitDir",boost::program_options::value<std::string>(&submitDir)->default_value("submitDir"),"name of output directory")
		("InDS,i",boost::program_options::value<std::string>(&InDS)->default_value(""),"InDS for grid job")
		("OutDS,o",boost::program_options::value<std::string>(&OutDS)->default_value(""),"OutDS for grid job")
		("Nevents,N",boost::program_options::value<int>(&Nevents)->default_value(-1),"Max number of events")
		("isGridJob",boost::program_options::value<bool>(&isGridJob)->default_value(false),"is it grid job?")
		("isTruth",boost::program_options::value<bool>(&istruth)->default_value(false),"Perform true analysis?")
		("isReco",boost::program_options::value<bool>(&isreco)->default_value(false),"Perform reco analysis?")
		("nFilesPerJob",boost::program_options::value<int>(&nFilesPerJob)->default_value(1),"Number of files per grid job")
		("reco_jet_collection",boost::program_options::value<std::string>(&reco_jet_collection)->default_value("HGamAntiKt4EMTopoJets"),"Jet collection")
		("ispTcut",boost::program_options::value<bool>(&ispTcut)->default_value(false),"Apply tight pT cuts on photons?")
		("isSignal",boost::program_options::value<bool>(&signal)->default_value(false),"Is it a signal sample?")
		("phPtcut",boost::program_options::value<float>(&phPtcut)->default_value(25),"Photon pT cut")
		("doCutFlow",boost::program_options::value<bool>(&doCutFlow)->default_value(false),"Do you want me to do the cutflow?")
		("doOptimization",boost::program_options::value<bool>(&optimization)->default_value(false),"Do you want to optimize the c tagging working point?")
		("jetpT",boost::program_options::value<float>(&jetpT)->default_value(25),"Jet pT cut")
		("jeteta",boost::program_options::value<float>(&jeteta)->default_value(2.5),"Jet #eta cut")
		("bfrac",boost::program_options::value<float>(&frac)->default_value(0.18),"The fraction of b")
		("DL1cut",boost::program_options::value<float>(&cut)->default_value(0.61456),"DL1 cut for c tagging")
		("btagger",boost::program_options::value<std::string>(&b_tagger)->default_value("MV2c10_FixedCutBeff_60"),"B tagger name")
		("dR_cut",boost::program_options::value<float>(&dR_cut)->default_value(0.5),"dR cut between a jet and photons")
		;


	//combine options types for parsing
	//all options may be specified on command line
	boost::program_options::options_description cmdline_options;
	cmdline_options.add(cmd_only_options).add(main_options);

	//all options except command line only may be specified in config file
	boost::program_options::options_description config_options;
	config_options.add(main_options);

	boost::program_options::variables_map vm;

	//first parse command line
	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdline_options), vm);
		boost::program_options::notify(vm);
	}
	catch(std::exception& e)
	{
		std::cerr << "Bad command line argument" << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	}

	//if config was specified, also parse config file
	if(vm.count("config"))
	{
		ifstream config_stream(config_file_name.c_str());
		try
		{
			boost::program_options::store(boost::program_options::parse_config_file(config_stream,cmdline_options), vm);
			boost::program_options::notify(vm);
		}
		catch(std::exception& e)
		{
			std::cerr << "Bad config file argument" << std::endl;
			std::cerr << e.what() << std::endl;
			return 1;
		}
	}

	//Need to be in MC mode for performance
	cout << endl << endl << "*********Configuration**********" << endl;
	if (Nevents>0) cout << "Number of events to be processed: " << Nevents <<endl;
	cout << "Input data set name for Grid: " << InDS <<endl;
	cout << "Output data set name for Grid: " << OutDS <<endl;
	cout << "Output file name in local: " << output_file_name <<endl;
	cout << "********************************" << endl << endl << endl;
	if(istruth){cout << "The NTupleAnalyser is configured to run on truth container" <<endl;}
	cout << "We will use Jet collection: " << reco_jet_collection << endl;
	if (doCutFlow) cout << "Initialising the cut-flow histogram" << endl;
	if (ispTcut) cout << "We will do for Tight photons a pT cut" << "which is: "<< phPtcut << endl;
	if (doCutFlow) cout << "Jet pT cut: " << jetpT << " Jet eta cut: " << jeteta << endl;
	if (signal) cout << "Applying the signal weights to MC"  << endl;
	cout << "dR cut between jet and photons is " << dR_cut << endl;



	// Set up the job for xAOD access:
	xAOD::Init().ignore();

	// Construct the samples to run on:
	SH::SampleHandler sh;

	// Get input file (! be careful about path -- not to include last folder !)
	if (!isGridJob){
		SH::ScanDir().filePattern("*").scan(sh, inputDir);
	}
	else {
		SH::scanDQ2 (sh, InDS.c_str());
		sh.setMetaString( "nc_grid_filter", "*");
	}
	sh.setMetaString( "nc_tree", "CollectionTree" );

	// Print what we found:
	sh.print();

	// Create an EventLoop job:
	cout << "Creating EventLoop job" << endl;
	EL::Job job;

	//Set outputFile
	EL::OutputStream output(output_file_name.c_str());
	job.outputAdd (output);
	EL::NTupleSvc *ntuple = new EL::NTupleSvc(output_file_name.c_str());
	job.algsAdd (ntuple);

	job.sampleHandler( sh );

	cout << "Seting maximum events to " << Nevents << endl;
	job.options()->setDouble (EL::Job::optMaxEvents, Nevents);

	// To automatically delete submitDir
	job.options()->setDouble(EL::Job::optRemoveSubmitDir, 1);


	// Add our analysis to the job:
	cout << "Add our analysis to the job" << endl;

	HcAnalyse* alg = new HcAnalyse();

    alg->_reco_jet_collection=reco_jet_collection.c_str();
    alg->_ispTcut = ispTcut;
    alg->_signal = signal;
    alg->_phpTcut = phPtcut;
    alg->_do_Cutflow = doCutFlow;
    alg->_jetpT = jetpT;
    alg->_jeteta = jeteta;
    alg->_frac = frac;
    alg->_cut = cut;
    alg->_Dooptimization = optimization;
    alg->_submitDir = submitDir;
    alg->_istruth = istruth;
    alg->_isreco = isreco;
    alg->_dR_cut = dR_cut;


	job.algsAdd( alg );
	alg->_outputName = output_file_name.c_str(); // give the name of the output to our algorithm
	// Run the job using the local/direct driver:
	cout << "Run the job" << endl;
	//Split level protection
    job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);

	if(!isGridJob)
	{
		 EL::DirectDriver driver;
		 driver.submit( job, submitDir );
	}
	else
	{
		EL::PrunDriver driver;
		driver.options()->setString("nc_outputSampleName",OutDS.c_str());
		if (nFilesPerJob > 0) driver.options()->setDouble("nc_nFilesPerJob", nFilesPerJob);
		driver.options()->setString("nc_cmtConfig", "x86_64-slc6-gcc62-opt");
		driver.options()->setDouble(EL::Job::optGridMergeOutput, 1); //run merging jobs for all samples before downloading (recommended)
		driver.submitOnly( job, submitDir );

	}
	cout << "We are done!" << endl;

	std::cout << "done looping" << std::endl;
}
