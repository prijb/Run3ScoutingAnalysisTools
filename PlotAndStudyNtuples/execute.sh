echo "Compiling the C++ files: "
g++ -I /home/arsahasransu/Software/madgraph/MG5_aMC_v2_7_0/HEPTools/boost/include/ robustanalyzermain.C data_robustanalyzer.C `root-config --cflags --glibs` -o data_robustanalyzer.out

echo "Compiling successful. Begin execution."

./data_robustanalyzer.out 0 &
proc0=$!
./data_robustanalyzer.out 1 &
proc1=$!
./data_robustanalyzer.out 2 &
proc2=$!
./data_robustanalyzer.out 3 &
proc3=$!
./data_robustanalyzer.out 4 &
proc4=$!
./data_robustanalyzer.out 5 &
proc5=$!
./data_robustanalyzer.out 6 &
proc6=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm data_robustanalyzer.out

hadd -f hists_DoubleElectronGunPt1To300.root hists_DoubleElectronGunPt1To300_?.root
hadd -f hists_DYToLLM4To50.root hists_DYToLLM4To50_?.root
hadd -f hists_ZeroBias2018D.root hists_ZeroBias2018D_?.root
hadd -f hists_Ephemeral1HLTPhysics2018D.root hists_Ephemeral1HLTPhysics2018D_?.root

rm hists_DoubleElectronGunPt1To300_?.root
rm hists_DYToLLM4To50_?.root
rm hists_ZeroBias2018D_?.root
rm hists_Ephemeral1HLTPhysics2018D_?.root
