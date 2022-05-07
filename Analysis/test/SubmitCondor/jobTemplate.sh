cd /user/asahasra/CMSSW_12_3_2/src/
eval `scram runtime -sh`

cd ./Run3ScoutingAnalysisTools/Analysis/test

ls
echo "=========== BEGINNING TO RUN THE SCRIPT ================"

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DoubleElectron_Pt-1To300-gun/ScoutingSkimSmall220504_DoubleElectronGunRun3Summer21_asahasra/220504_100553/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=120X_mcRun3_2021_realistic_v6 output=DoubleElectronGunPt1To300_ScoutingSkimSmall220504_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DoubleElectron_Pt-1To300-gun/ScoutingSkimSmall220504_DoubleElectronGunRun3Summer21_asahasra/220504_100553/0000/outputScoutingPF_$num.root
done

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkimSmall220504_QCDPt20To30EMEnrichedRun3Summer21_asahasra/220504_115849/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=120X_mcRun3_2021_realistic_v6 output=QCDPt20To30EmEnriched_ScoutingSkim220504_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkimSmall220504_QCDPt20To30EMEnrichedRun3Summer21_asahasra/220504_115849/0000/outputScoutingPF_$num.root
done

echo "================= END OF SCRIPT RUN ===================="
ls
