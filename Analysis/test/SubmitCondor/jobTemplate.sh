cd /user/asahasra/ProcessScoutingJPsi/CMSSW_13_0_0_pre3/src/
eval `scram runtime -sh`

cd ./Run3ScoutingAnalysisTools/Analysis/test

ls
echo "=========== BEGINNING TO RUN THE SCRIPT ================"

for num in $(ls /pnfs/iihe/cms/store/user/asahasra/JPsiToEE_Pt2To30_13p6TeV_TuneCP5_pythia8/PrivateMC_CMSSW1300pre3_230207ScoutingReRun/230207_120037/0000/ | grep -o -E '[0-9]+')
do
    # Uncomment if file contains 2022 in its name
    #if [[ $num -eq 2022 ]]; then
    #   continue
    #fi
    # Check if $num ends with $1
    if [[ $(($num%10)) -ne $1 ]]; then
	continue
    fi
    cmsRun ScoutingNanoAOD_cfg.py isMC=True GlobalTagMC=124X_mcRun3_2022_realistic_forTSG_menu1p4_v1 output=Scouting_JPsiToEE_nTuple_220207Full_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/JPsiToEE_Pt2To30_13p6TeV_TuneCP5_pythia8/PrivateMC_CMSSW1300pre3_230207ScoutingReRun/230207_120037/0000/outputScoutingPF_$num.root
done

echo "================= END OF SCRIPT RUN ===================="
ls
