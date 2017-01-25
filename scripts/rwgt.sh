############################################
#                                          #
# Script to reweight the tarred event files #
#                                          #
############################################

RWGT=""

tarred=1

# Rename all scripts (just in case they aren't in ascending order)

if [ -e events ] ; then
    rm events
fi

ls pwgevents*.gz >> events

numFiles=`ls pwgevents*.gz | wc -l` # this is how many files we have to rwgt
numProcessors=$((numFiles/4)) # how many processors we need in total (will split this up into each script doing loops rather than submitting large numbers of scripts)
numScripts=30 # this is how many scripts we want to runnumScripts=$((numFiles/4)) # this is how many scripts we need (each script runs on 4 processors)
bigLoops=$((numProcessors/numScripts))
if [ ! $((numProcessors % numScripts)) -eq 0 ] ; then
    bigLoops=$((bigLoops+1)) # add an extra one to pick up files which don't divide evenly
fi

for i in `seq 1 $numScripts` ; do

    file=`basename $PWD`-rwgt-$i.pbs
    cp ../scripts/rwgt.pbs $file
    sed -i "s/nScripts=.*/nScripts=$numScripts/g" $file
    sed -i "s/scriptNumber=.*/scriptNumber=$i/g" $file
    sed -i "s/maxBigLoops=.*/maxBigLoops=$bigLoops/g" $file


    if [ $i -eq 1 ] ; then
        RWGT=$(qsub $file)
    else
        RWGT=$RWGT:$(qsub $file)
    fi

done
