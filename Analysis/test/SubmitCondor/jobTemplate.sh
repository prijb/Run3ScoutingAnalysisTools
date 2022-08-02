#!/bin/sh

cd $TMPDIR
mkdir Job_$1_$2
ls
echo "============ Created the master run folder ============"

cd Job_$1_$2
cd /afs/cern.ch/work/a/asahasra/public/ScoutingStudy/Process13p6data/CMSSW_12_4_4/src/
eval `scramv1 runtime -sh`
cd -

cp /afs/cern.ch/work/a/asahasra/public/ScoutingStudy/Process13p6data/CMSSW_12_4_4/src/Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py ./

ls
echo "============ Copied the python config file ============"
echo $3
#cp $3 ./
echo "============ Listing the directory contents pre cmsRun ============"
ls

echo "============ Running ============"

cmsRun ScoutingNanoAOD_cfg.py isMC=False GlobalTagData=124X_dataRun3_HLT_v4 output=dataNtuple_$1_$2.root inputFile=file:$3

echo "============ Listing the directory contents post cmsRun ============"
ls

cp dataNtuple_$1_$2.root /eos/user/a/asahasra/Scouting13p6TeVEGammaNtuples/Run2022C_Skim220801
echo "============ All processes completed. Exiting now.  ============"
