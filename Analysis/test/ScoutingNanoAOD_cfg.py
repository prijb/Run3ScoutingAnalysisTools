import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'isMC',
    True,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample is simulation or data'
)

params.register(
    'useWeights',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'filterTrigger',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to fire a trigger used in the analysis'
)

params.register(
    'filterMuons',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to contain at least two muons'
)

params.register(
    'reducedInfo',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to store just the reduced information'
)

params.register(
    'trigProcess',
    'HLT',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagData', 
    '112X_mcRun3_2021_realistic_v16', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagMC',
    '112X_mcRun3_2021_realistic_v16',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'xsec',
    0.001,
    VarParsing.multiplicity.singleton,VarParsing.varType.float,
    'Cross-section for a Monte Carlo Sample'
)

params.register(
    'output', 
    'scoutingNTuple.root', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Name of the output root file'
)

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120001/003f1472-9e77-47d9-9ea4-5dec6dbc50ad.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/52da57ee-fdb4-4268-b8ff-c4b23da87d2b.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/64c37853-fe4e-4440-a23b-6ae4fcdd23f0.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/87e3862d-c57f-43da-9763-e4d6fcb6b93b.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/e5dd5c81-1486-4538-8a3e-5e887781744f.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/f0da900d-a175-4e6a-8540-91b6e47f945b.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/9d650050-bd6d-42c4-8dab-2e4ed4bb8b86.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/402501fc-9cd1-4613-aef3-1cda21eee13e.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/799a2712-2996-4e76-832e-ea2fdf8f265c.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/81a15d64-5617-49f0-915b-e5f4ab692ceb.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/28f5a6b8-149a-4cca-93c9-89e6787fdeb9.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/907b8e71-05e4-484d-9848-fc309a14f30a.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/9591747e-abf4-4743-8bb0-d4ed81eb6b7b.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/7e581e35-7cc5-4bd9-922b-b72c93f03b88.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/95696d78-a123-4c5f-a0e4-23afb4902c88.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/262cb0c4-1e74-4591-9220-c7560b89a8b3.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/d7c102d3-50d0-41d6-a03d-b9509cec9e93.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/425b166b-4baf-4d55-aaf5-0c8eabfda3c9.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/8ec3011e-6168-463d-9bd3-5550a7bc1d4b.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/703c9678-9b11-4a92-b803-5538021afabe.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/e0af59e0-3d24-410f-bf75-c33fdcf8806b.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/d5e5e7d7-2ae8-42ea-871b-67089da81129.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280009/42a645fb-b838-4388-9a3c-a90deef4de32.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120011/0ffb0b1c-51eb-4bc9-9e62-ac06bcaa25d0.root',
                                'root://xrootd-cms.infn.it//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-RECO/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/120011/c0fc4976-fd04-4382-9246-e4ecb17a705e.root'
                            ])
                        )

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.isMC :
    process.GlobalTag.globaltag = params.GlobalTagMC
else :
    process.GlobalTag.globaltag = params.GlobalTagData

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(params.output)
                               )
#from DarkPhotonAnalysis.DimuonAnalysis2018.TriggerPaths_cfi import getL1Conf
L1Info = ['L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu_12_5','L1_DoubleMu_15_7','L1_TripleMu_5_3_3','L1_TripleMu_5_5_3','L1_QuadMu0','L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4','L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18','L1_DoubleMu4_SQ_OS_dR_Max1p2','L1_SingleMu22','L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4','L1_DoubleMu4p5_SQ_OS_dR_Max1p2','L1_DoubleMu4p5_SQ_OS','L1_DoubleMu0er1p5_SQ_dR_Max1p4','L1_DoubleMu0er2p0_SQ_dR_Max1p4','L1_DoubleMu0_SQ']
# Make tree
process.mmtree = cms.EDAnalyzer('EGammaOnly_ScoutingNanoAOD',
                                triggerresults   = cms.InputTag("TriggerResults", "", params.trigProcess),
                                doL1 = cms.bool(False),
                                triggerConfiguration = cms.PSet(
                                    hltResults            = cms.InputTag('TriggerResults','','HLT'),
                                    l1tResults            = cms.InputTag(''),
                                    daqPartitions         = cms.uint32(1),
                                    l1tIgnoreMaskAndPrescale = cms.bool(False),
                                    throw                 = cms.bool(False)
                                ),
                                ReadPrescalesFromFile = cms.bool( False ),
                                AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1Seeds           = cms.vstring(L1Info),
                                beamspot        = cms.InputTag("offlineBeamSpot"),
                                electrons        = cms.InputTag("hltScoutingEgammaPacker"),
                                photons          = cms.InputTag("hltScoutingEgammaPacker"),
                                gens = cms.InputTag("genParticles"),
                                
                            )

process.p = cms.Path(                  process.mmtree)
