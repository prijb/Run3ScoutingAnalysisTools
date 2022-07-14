#!/bin/sh

cd $TMPDIR
mkdir Job_$1_$2
ls
echo "============ Created the master run folder ============"

cd Job_$1_$2
cd /afs/cern.ch/work/a/asahasra/private/ScoutingPFRun3_RAW_v1_000_353_706/CMSSW_12_4_0/src
eval `scramv1 runtime -sh`
cd -

cp /afs/cern.ch/work/a/asahasra/private/ScoutingPFRun3_RAW_v1_000_353_706/CMSSW_12_4_0/src/Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py ./

ls
echo "============ Copied the python config file ============"
echo $3
#cp $3 ./
echo "============ Listing the directory contents pre cmsRun ============"
ls

echo "============ Running ============"

cmsRun ScoutingNanoAOD_cfg.py isMC=False GlobalTagData=123X_dataRun3_HLT_v7 output=dataNtuple_$1_$2.root inputFile=file:$3

echo "============ Listing the directory contents post cmsRun ============"
ls

cp dataNtuple_$1_$2.root /afs/cern.ch/work/a/asahasra/private/ScoutingPFRun3_RAW_v1_000_353_706/CMSSW_12_4_0/src/Run3ScoutingAnalysisTools/Analysis/test/SubmitCondor/
echo "============ All processes completed. Exiting now.  ============"
