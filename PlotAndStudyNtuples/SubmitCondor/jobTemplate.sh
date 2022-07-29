#!/bin/sh

cd $TMPDIR
mkdir Job_$1_$2
echo "============ Created the master run folder ============"

cd Job_$1_$2
cd /afs/cern.ch/work/a/asahasra/public/ScoutingStudy/Process13p6data/CMSSW_12_4_0/src/
eval `scramv1 runtime -sh`
cd $TMPDIR/Job_$1_$2
pwd

cp $CMSSW_BASE/src/Run3ScoutingAnalysisTools/PlotAndStudyNtuples/robustanalyzer* ./

echo "============ Copied the cpp run files file. Listing directory contents ============"
ls
echo "======================= End listing the directory contents ========================"
echo ""

echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out $2 49

echo ""
echo "Run over. Copy the output files."

cp hists_data_*.root $CMSSW_BASE/src/Run3ScoutingAnalysisTools/PlotAndStudyNtuples/hists_data_$1_$2.root
