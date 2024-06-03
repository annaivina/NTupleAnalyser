# The name of the package:
atlas_subdir (NTupAnalyser)


# Component(s) in the package:
atlas_add_root_dictionary (NTupAnalyserLib NTupAnalyserCintDict
   ROOT_HEADERS NTupAnalyser/*.h Root/LinkDef.h
   EXTERNAL_PACKAGES ROOT )

# add the library itself:
atlas_add_library (NTupAnalyserLib
   NTupAnalyser/*.h  Root/*.cxx ${NTupAnalyserCintDict}
   PUBLIC_HEADERS NTupAnalyser
   INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS}
   LINK_LIBRARIES ${EIGEN_LIBRARIES}
   AsgTools AthContainers AthLinks GeoPrimitives EventLoopGrid EventLoop EventLoopAlgs PathResolver GoodRunsListsLib TrigDecisionToolLib xAODCaloEvent xAODEventInfo xAODTracking xAODTruth xAODHIEvent xAODTrigMinBias xAODJet xAODMuon xAODForward TrigConfxAODLib SampleHandler InDetTrackSelectionToolLib ElectronPhotonSelectorToolsLib ElectronPhotonFourMomentumCorrectionLib MuonSelectorToolsLib MuonMomentumCorrectionsLib JetCalibToolsLib xAODEgamma
   PRIVATE_LINK_LIBRARIES )


# The dependencies of the package:
atlas_depends_on_subdirs(
   PUBLIC
   Control/AthContainers
   Control/AthToolSupport/AsgTools
   Control/xAODRootAccessInterfaces
   Control/xAODRootAccess
   Control/CxxUtils
   Event/EventPrimitives
   Event/xAOD/xAODBase
   Event/xAOD/xAODEgamma
   Event/xAOD/xAODJet
   Event/xAOD/xAODMuon
   Event/xAOD/xAODTau
   Event/xAOD/xAODCore
   Event/xAOD/xAODMissingET
   Event/xAOD/xAODEventInfo
   Event/xAOD/xAODTruth
   Event/xAOD/xAODForward
   Event/xAOD/xAODCaloEvent
   Event/xAOD/xAODTracking
   Event/xAOD/xAODTrigMinBias
   Event/xAOD/xAODTrigger
   Event/xAOD/xAODHIEvent
   Event/xAOD/xAODRootAccess
   DataQuality/GoodRunsLists/GoodRunsListSelectionTool
   InnerDetector/InDetRecTools/InDetTrackSelectionTool
   PhysicsAnalysis/AnalysisCommon/AssociationUtils
   PhysicsAnalysis/AnalysisCommon/IsolationSelection
   PhysicsAnalysis/AnalysisCommon/PATInterfaces
   PhysicsAnalysis/D3PDTools/RootCoreUtils
   PhysicsAnalysis/D3PDTools/SampleHandler
   PhysicsAnalysis/D3PDTools/EventLoop
   PhysicsAnalysis/AnalysisCommon/PileupReweighting
   PhysicsAnalysis/ElectronPhotonID/ElectronEfficiencyCorrection
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonShowerShapeFudgeTool
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonSelectorTools
   PhysicsAnalysis/ElectronPhotonID/IsolationCorrections
   PhysicsAnalysis/ElectronPhotonID/PhotonEfficiencyCorrection
   PhysicsAnalysis/JetMissingEtID/JetSelectorTools
   PhysicsAnalysis/JetTagging/JetTagPerformanceCalibration/xAODBTaggingEfficiency
   PhysicsAnalysis/MuonID/MuonsIDAnalysis/MuonEfficiencyCorrections
   PhysicsAnalysis/MuonID/MuonsIDAnalysis/MuonMomentumCorrections
   PhysicsAnalysis/MuonID/MuonSelectorTools
   PhysicsAnalysis/TauID/TauAnalysisTools
   Reconstruction/Jet/JetCalibTools
   Reconstruction/Jet/JetJvtEfficiency
   Reconstruction/Jet/JetMomentTools
   Reconstruction/Jet/JetResolution
   Reconstruction/Jet/JetUncertainties
   Reconstruction/MET/METInterface
   Reconstruction/MET/METUtilities
   Tools/PathResolver
   Trigger/TrigAnalysis/TrigDecisionTool
   Trigger/TrigAnalysis/TriggerMatchingTool
   Trigger/TrigConfiguration/TrigConfxAOD )

# Extra dependencies, based on the build environment:
#set( extra_deps )
#if( XAOD_STANDALONE )
#    set( extra_deps PhysicsAnalysis/D3PDTools/EventLoop
#    PRIVATE Control/xAODRootAccessInterfaces Control/xAODRootAccess )
#else()
#    set( extra_deps PRIVATE Control/AthenaBaseComps Control/AthAnalysisBaseComps
#    Control/StoreGate GaudiKernel )
#endif()

atlas_install_data( share/*.xml
   DESTINATION share )

# External(s) used by the package:
find_package( Boost COMPONENTS program_options )
find_package( ROOT COMPONENTS Core Hist RIO Gpad Tree MathCore Graf )
find_package( Eigen )


atlas_add_executable( run_analyser
   util/run_analyser.cxx 
   LINK_LIBRARIES NTupAnalyserLib )

atlas_add_executable( run_ntuple
   util/run_ntuple.cxx 
   LINK_LIBRARIES NTupAnalyserLib )


