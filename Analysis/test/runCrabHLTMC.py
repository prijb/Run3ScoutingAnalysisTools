from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3Scouting_STHDM_M200dM20ctau3cm_Run3Winter2021MC_crabRunSkim2'
config.General.workArea = 'crab_Run3Scouting_STHDM_M200dM20ctau3cm_Run3Winter2021MC'
config.General.transferOutputs = True
config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ScoutingNanoAOD_cfg.py'

config.Data.inputDataset = '/SingletTripletHDMToDisplacedL_M200deltaM20ctau3cm_TuneCP5_14TeV-madgraph-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_rndm_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.publication = False
config.Data.outputDatasetTag = 'STHDM_M200dM20ctau3cm_Run3Winter2021MC_crabRunSkim2_asahasra'
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T2_BE_IIHE'
