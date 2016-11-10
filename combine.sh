cd combined-output
for i in `seq 1 1` ; do
    mkdir combined_W$i
    cp ../LHEF-output/pwgLHEF*W$i.top combined_W$i
    cp ../PYTHIA8-output/pwg-*POWHEG+PYTHIA8-output-W$i.top combined_W$i
    cp ../mergedata combined_W$i
    cd combined_W$i
    ./mergedata 1 pwgLHEF*top
    rm pwgLHEF*top
    mv fort.12 LHEF-W$i.top
    ./mergedata 1 pwg-*POWHEG+PYTHIA8-output-W$i.top
    rm pwg*top
    mv fort.12 P8-W$i.top
    cd ../
done
