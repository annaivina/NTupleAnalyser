#include <EventLoop/DirectDriver.h>
#include <EventLoop/Job.h>

#include <TSystem.h>
#include <SampleHandler/ScanDir.h>
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/OutputStream.h"
#include "EventLoopAlgs/NTupleSvc.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <iostream>
#include "xAODRootAccess/Init.h"


#include "NTupAnalyser/HcAnalyse.h"

#include <iostream>



// to submit to grid next time:
// run_overlay_treemaker 6 1

int main(int argc, char *argv[])
{
	char submitDir[100]="submitDir";
	int vers=0;
	int isGridJob=0;

	if( argc > 1 )
	{
		vers=atoi(argv[1]);
		sprintf(submitDir,"submitDir_v%d", vers);
	}

	if( argc > 2 )
	{
		isGridJob=atoi(argv[2]);
	}

	xAOD::Init().ignore();
	SH::SampleHandler sh;
	char filePath[300];

	// Set up the job for xAOD access:
	xAOD::Init().ignore();

	if(!isGridJob)
	{
		switch(vers)
		{
			case 1: sprintf(filePath, ""); break;

			default: std::cout<<"fail "<<vers<<std::endl; return -1;
		}

		const char* inputFilePath = gSystem->ExpandPathName(filePath);
		SH::ScanDir().filePattern("*AOD*").scan(sh, inputFilePath);
	}
	else
	{
		SH::scanDQ2(sh,"data17_hi.00338037.physics_MinBiasOverlay.merge.AOD.f900_m1912");
		sh.setMetaString( "nc_grid_filter", "*AOD*");
	}
	sh.setMetaString( "nc_tree", "CollectionTree" );
	sh.print();

	EL::Job job;
	job.sampleHandler( sh );
	job.options()->setDouble (EL::Job::optMaxEvents, -1);
	job.options()->setDouble(EL::Job::optRemoveSubmitDir, 1);
	//job.options()->setDouble(EL::Job::optGridNGBPerJob, 2000);
	job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);
	job.options()->setDouble(EL::Job::optGridNFilesPerJob, 15);


	HcAnalyse* alg = new HcAnalyse();
	job.algsAdd(alg);



	if(!isGridJob)
	{
		 EL::DirectDriver driver;
		 driver.submit( job, submitDir );
	}
	else
	{
		EL::PrunDriver driver;
		driver.options()->setString("nc_outputSampleName",Form("user.pbalek.overlay_XeXe.338037.v%d",vers));
		//driver.options()->setDouble("nc_nFilesPerJob", nFilesPerJob);
		driver.options()->setString("nc_cmtConfig", "x86_64-slc6-gcc62-opt");
		driver.options()->setDouble(EL::Job::optGridMergeOutput, 1);

		driver.submitOnly( job, submitDir );
	}


}

