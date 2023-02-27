#!/bin/sh

THISP=$(pwd)
mkdir Job_$1_$2
ls
echo "============ Created the master run folder ============"

cd /user/asahasra/ProcessScoutingJPsi/CMSSW_13_0_0_pre3/src/
eval `scram runtime -sh`

cd $THISP
cd Job_$1_$2

cp /user/asahasra/ProcessScoutingJPsi/CMSSW_13_0_0_pre3/src/Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py ./

ls
echo "============ Copied the python config file ============"
echo $3
#cp $3 ./
echo "============ Listing the directory contents pre cmsRun ============"
ls
echo "=========== BEGINNING TO RUN THE SCRIPT ================"

cmsRun ScoutingNanoAOD_cfg.py isMC=False GlobalTagData=126X_dataRun3_HLT_v1 output=ScoutingNTuple_Eph02022G_230224trial_$1_$2.root inputFile=file:$3

echo "================= END OF SCRIPT RUN ===================="
ls

cp ScoutingNTuple_Eph02022G_230224trial_$1_$2.root /user/asahasra/ProcessScoutingJPsi/CMSSW_13_0_0_pre3/src/Run3ScoutingAnalysisTools/Analysis/test/
echo "============ All processes completed. Exiting now.  ============"
