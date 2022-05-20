cd /user/asahasra/CMSSW_12_3_2/src/
eval `scram runtime -sh`

cd ./Run3ScoutingAnalysisTools/Analysis/test

ls
echo "=========== BEGINNING TO RUN THE SCRIPT ================"

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DoubleElectron_Pt-1To300-gun/ScoutingSkimSmall220518_2_DoubleElectronGunRun3Summer21_asahasra/220518_132349/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True GlobalTagMC=123X_mcRun3_2021_realistic_v4 output=DoubleElectronGunPt1To300_ScoutingSkimSmall220518_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DoubleElectron_Pt-1To300-gun/ScoutingSkimSmall220518_2_DoubleElectronGunRun3Summer21_asahasra/220518_132349/0000/outputScoutingPF_$num.root
done
for num in $(ls /pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkimSmall220518_QCDPt20To30EMEnrichedRun3Summer21_asahasra/220518_131330/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True GlobalTagMC=123X_mcRun3_2021_realistic_v4 output=QCDPt20To30EMEnriched_ScoutingSkimSmall220518_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkimSmall220518_QCDPt20To30EMEnrichedRun3Summer21_asahasra/220518_131330/0000/outputScoutingPF_$num.root
done
for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkimSmall220518_DYToLLM4To50Run3Summer21_asahasra/220519_072112/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True GlobalTagMC=123X_mcRun3_2021_realistic_v4 output=DYToLLM4To50_ScoutingSkimSmall220518_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkimSmall220518_DYToLLM4To50Run3Summer21_asahasra/220519_072112/0000/outputScoutingPF_$num.root
done
for num in $(ls /pnfs/iihe/cms/store/user/asahasra/ZeroBias/ScoutingSkimSmall220518_ZeroBias2018D_asahasra/220518_130504/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=False GlobalTagData=123X_dataRun3_HLT_frozen_v2 output=ZeroBias2018D_ScoutingSkimSmall220518_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/ZeroBias/ScoutingSkimSmall220518_ZeroBias2018D_asahasra/220518_130504/0000/outputScoutingPF_$num.root
done
echo "================= END OF SCRIPT RUN ===================="
ls
