############################################
#                                          #
# Script to analyse the tarred event files #
#                                          #
############################################

do_lhef=1
do_py8=1

# Rename all scripts (just in case they aren't in ascending order)

if [ -e events ] ; then
    rm events
fi

ls pwgevents*.gz >> events

numFiles=`ls pwgevents*.gz | wc -l` # this is how many files we have to analyse
numScripts=$((numFiles/4)) # this is how many scripts we need (each script runs on 4 processors)

for i in `seq 1 $numScripts` ; do

    file=`basename $PWD`-analyse-$i.pbs
    cp ../scripts/analyse.pbs $file
    sed -i "s/nScripts=.*/nScripts=$numScripts/g" $file
    sed -i "s/scriptNumber=.*/scriptNumber=$i/g" $file

    if [ $do_lhef -eq 1 ] ; then
        sed -i "s/do_lhef=.*/do_lhef=1/g" $file
    else
        sed -i "s/do_lhef=.*/do_lhef=0/g" $file
    fi

    if [ $do_py8 -eq 1 ] ; then
        sed -i "s/do_py8=.*/do_py8=1/g" $file
    else
        sed -i "s/do_py8=.*/do_py8=0/g" $file
    fi

    qsub $file

done