import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'inputFile', 
    'HLT2022_HLT_OR_GENSIMDIGIRAW.root', 
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Name of the input root file'
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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(params.inputFile),
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
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(params.output)
                               )

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
                            )

#process.p = cms.Path( process.gtStage2Digis*process.mmtree )
process.p = cms.Path( process.mmtree )
