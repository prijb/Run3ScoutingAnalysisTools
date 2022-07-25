echo "Compiling the C++ files: "
g++ -I /home/arsahasransu/Software/madgraph/MG5_aMC_v2_7_0/HEPTools/boost/include/ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 &
proc0=$!
./robustanalyzer.out 1 &
proc1=$!
./robustanalyzer.out 2 &
proc2=$!
./robustanalyzer.out 3 &
proc3=$!
./robustanalyzer.out 4 &
proc4=$!
./robustanalyzer.out 5 &
proc5=$!
./robustanalyzer.out 6 &
proc6=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm robustanalyzer.out

hadd -f hists_data.root hists_data_?.root
#hadd -f hists_DoubleElectronGunPt1To300.root hists_DoubleElectronGunPt1To300_?.root
#hadd -f hists_DoubleElectronGunPt1To300Old.root hists_DoubleElectronGunPt1To300Old_?.root
#hadd -f hists_QCDPt20To30EmEnriched.root hists_QCDPt20To30EmEnriched_?.root
#hadd -f hists_QCDPt30To50EmEnriched.root hists_QCDPt30To50EmEnriched_?.root

rm hists_data_?.root
#rm hists_DoubleElectronGunPt1To300_?.root
#rm hists_DoubleElectronGunPt1To300Old_?.root
#rm hists_QCDPt20To30EmEnriched_?.root
#rm hists_QCDPt30To50EmEnriched_?.root
