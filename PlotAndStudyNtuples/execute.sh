echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 6 &
proc0=$!
./robustanalyzer.out 1 6 &
proc1=$!
./robustanalyzer.out 2 6 &
proc2=$!
./robustanalyzer.out 3 6 &
proc3=$!
./robustanalyzer.out 4 6 &
proc4=$!
./robustanalyzer.out 5 6 &
proc5=$!
./robustanalyzer.out 6 6 &
proc6=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" ]
#while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" ]
do
    echo "====Still executing===="
    sleep 5
done

echo "Run over. Clean up and combine files."

rm robustanalyzer.out

hadd -f hists_Eph02022G.root hists_Eph02022G_?.root
#hadd -f hists_JPsi.root hists_JPsi_?.root

rm hists_Eph02022G_?.root
#rm hists_JPsi_?.root
