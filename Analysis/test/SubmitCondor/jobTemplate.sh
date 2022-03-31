cd /user/asahasra/CMSSW_12_1_0_pre4/src/
eval `scram runtime -sh`

cd ./Run3ScoutingAnalysisTools/Analysis/test

ls
echo "=========== BEGINNING TO RUN THE SCRIPT ================"

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DYToLL_M-50_TuneCP5_14TeV-pythia8/ScoutingSkim220201_DYToLLM50Run3Summer21_asahasra/220201_162520/0000/ | grep -o -E '[0-9]+')
do
    if [[ $num -eq 2022 ]]; then
	continue
    fi
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=DYToLLM50_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DYToLL_M-50_TuneCP5_14TeV-pythia8/ScoutingSkim220201_DYToLLM50Run3Summer21_asahasra/220201_162520/0000/HLT2022_HLT_$num.root
done

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkim220127_DYToLLM4To5Run3Summer21_asahasra/220128_130307/0000/ | grep -o -E '[0-9]+')
do
    if [[ $num -eq 2022 ]]; then
	continue
    fi
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=DYToLLM4To50_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkim220127_DYToLLM4To5Run3Summer21_asahasra/220128_130307/0000/HLT2022_HLT_$num.root
done

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/ | grep -o -E '[0-9]+')
do
    if [[ $num -eq 2022 ]]; then
	continue
    fi
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=QCDPt20To30EmEnriched_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_$num.root
done

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt30To50EmEnrichedRun3Summer21_asahasra/220128_130341/0000/ | grep -o -E '[0-9]+')
do
    if [[ $num -eq 2022 ]]; then
	continue
    fi
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=QCDPt30To50EmEnriched_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt30To50EmEnrichedRun3Summer21_asahasra/220128_130341/0000/HLT2022_HLT_$num.root
done

hadd -f DYToLLM50_ScoutingSkim220127_$1.root DYToLLM50_ScoutingSkim220127_$1_*.root
hadd -f DYToLLM4To50_ScoutingSkim220127_$1.root DYToLLM4To50_ScoutingSkim220127_$1_*.root
hadd -f QCDPt20To30EmEnriched_ScoutingSkim220127_$1.root QCDPt20To30EmEnriched_ScoutingSkim220127_$1_*.root
hadd -f QCDPt30To50EmEnriched_ScoutingSkim220127_$1.root QCDPt30To50EmEnriched_ScoutingSkim220127_$1_*.root

rm DYToLLM50_ScoutingSkim220127_$1_*.root
rm DYToLLM4To50_ScoutingSkim220127_$1_*.root
rm QCDPt20To30EmEnriched_ScoutingSkim220127_$1_*.root
rm QCDPt30To50EmEnriched_ScoutingSkim220127_$1_*.root

echo "================= END OF SCRIPT RUN ===================="
ls
