import FWCore.ParameterSet.Config as cms

# Define the process
process = cms.Process("LL")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 10 ) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/360/390/00000/c4c4cf1c-ca37-4aa3-838a-78638a5d296e.root"),
                            secondaryFileNames = cms.untracked.vstring(
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/14843255-6295-4a01-a224-d3411f2dd7cd.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/4baea980-363a-4a93-81dc-8960bd8c5551.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/5b5073cd-6600-4081-a345-a83ef75bc64b.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/ae50cca3-45b9-49d0-a563-84e6dc6591e4.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/cccae9d4-33c8-4519-81bf-aecdebbac390.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/e116cb80-e579-4138-ac12-505b0a3986e9.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/f6c3952b-f1ea-4b58-9173-bf35e304029b.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/f9a71e47-6c0a-4710-9741-1857acc378c7.root",
                                "root://xrootd-cms.infn.it///store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/360/390/00000/fe3ada5c-a442-4725-9504-5cdfeac14e9d.root"
                            ),
                        )

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '124X_dataRun3_Prompt_v10'

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('scoutingNTuple.root')
                               )
L1Info = ['L1_DoubleMu_12_5', 'L1_DoubleMu_15_7', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18', 'L1_DoubleMu4_SQ_OS_dR_Max1p2', 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2', 'L1_HTT200er', 'L1_HTT255er', 'L1_HTT280er', 'L1_HTT320er', 'L1_HTT360er', 'L1_HTT400er', 'L1_HTT450er', 'L1_ETT2000', 'L1_SingleJet180', 'L1_SingleJet200', 'L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5', 'L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5', 'L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5', 'L1_SingleLooseIsoEG28er2p1', 'L1_SingleLooseIsoEG28er1p5', 'L1_SingleLooseIsoEG30er1p5', 'L1_SingleIsoEG28er2p1', 'L1_SingleIsoEG30er2p1', 'L1_SingleIsoEG32er2p1', 'L1_DoubleEG_LooseIso16_LooseIso12_er1p5', 'L1_DoubleEG_LooseIso18_LooseIso12_er1p5', 'L1_DoubleEG_LooseIso20_LooseIso12_er1p5', 'L1_DoubleEG_LooseIso22_LooseIso12_er1p5']

# Make tree
process.mmtree = cms.EDAnalyzer('EGammaOnly_ScoutingNanoAOD',
                                triggerhltres = cms.InputTag("TriggerResults::HLT"),
                                oflprimaryVtx = cms.InputTag("offlineSlimmedPrimaryVerticesWithBS::RECO"),
                                oflelectrons = cms.InputTag("slimmedElectrons::RECO"),
                                oflrho = cms.InputTag("fixedGridRhoAll::RECO"),
                                primaryVtx = cms.InputTag("hltScoutingPrimaryVertexPacker:primaryVtx:HLT"),
                                electrons = cms.InputTag("hltScoutingEgammaPacker::HLT"),
                                photons = cms.InputTag("hltScoutingEgammaPacker::HLT"),
                                rho = cms.InputTag("hltScoutingPFPacker:rho"),
                                AlgInputTag = cms.InputTag("gtStage2Digis"),
                                l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1Seeds = cms.vstring(L1Info),
                                doL1 = cms.bool(True)
                            )

process.p = cms.Path( process.mmtree )
