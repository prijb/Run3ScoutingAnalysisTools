from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3Scouting2022F_crabRunSkim230502'
config.General.workArea = 'crabRunSkim230502_Run3Scouting2022F'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ScoutingNanoAOD_forCrab_cfg.py'
config.JobType.inputFiles = ['runCrabScoutingMonitorData2022F.py']

config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-PromptReco-v1/MINIAOD'
config.Data.useParent = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.LumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
config.Data.publication = False
config.Data.outputDatasetTag = 'Scouting2022F_crabRunSkim230502_asahasra'
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T2_BE_IIHE'
