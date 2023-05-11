#!/bin/sh

cd $TMPDIR
mkdir Job_$1_$2
ls
echo "============ Created the master run folder ============"

cd Job_$1_$2
cd /user/asahasra/CommissioningPromptElectrons2023/CMSSW_13_0_5_patch1/src/
eval `scramv1 runtime -sh`
cd -

cp /user/asahasra/CommissioningPromptElectrons2023/CMSSW_13_0_5_patch1/src/Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py ./

ls
echo "============ Copied the python config file ============"
echo $3
#cp $3 ./
echo "============ Listing the directory contents pre cmsRun ============"
ls

echo "============ Running ============"

cmsRun ScoutingNanoAOD_cfg.py GlobalTagData=130X_dataRun3_HLT_v2 inputFile=file:$3 output=scoutingNTuple_$1_$2.root

echo "============ Listing the directory contents post cmsRun ============"
ls

cp scoutingNTuple_$1_$2.root /user/asahasra/CommissioningPromptElectrons2023/CMSSW_13_0_5_patch1/src/Run3ScoutingAnalysisTools/Analysis/test
echo "============ All processes completed. Exiting now.  ============"
