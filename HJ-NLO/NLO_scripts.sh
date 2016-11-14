startingSeed=1
numScripts=20

# Stage 1 - xgrid 1
XG1=""
for i in `seq 1 $numScripts` ; do
    cp ../scripts/st1-1.pbs `basename $PWD`-xg1-$i.pbs
    sed -i "s/nScripts=15/nScripts=$numScripts/g" `basename $PWD`-xg1-$i.pbs
    sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-xg1-$i.pbs
    if [ $i -eq 1 ] ; then
	XG1=$(qsub `basename $PWD`-xg1-$i.pbs)
    else
	XG1=$XG1:$(qsub `basename $PWD`-xg1-$i.pbs)
    fi
done
# Stage 1 - xgrid 2
XG2=""
for i in `seq 1 $numScripts` ; do
    cp ../scripts/st1-2.pbs `basename $PWD`-xg2-$i.pbs
    sed -i "s/nScripts=15/nScripts=$numScripts/g" `basename $PWD`-xg2-$i.pbs
    sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-xg2-$i.pbs
    if [ $i -eq 1 ] ; then
	XG2=$(qsub -W depend=afterany:$XG1 `basename $PWD`-xg2-$i.pbs)
    else
	XG2=$XG2:$(qsub -W depend=afterany:$XG1 `basename $PWD`-xg2-$i.pbs)
    fi
done

# Stage 2 - NLO and btilde upper bound
ST2=""
for i in `seq 1 $numScripts` ; do
    cp ../scripts/st2.pbs `basename $PWD`-st2-$i.pbs
    sed -i "s/nScripts=15/nScripts=$numScripts/g" `basename $PWD`-st2-$i.pbs
    sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-st2-$i.pbs
    if [ $i -eq 1 ] ; then
	ST2=$(qsub -W depend=afterany:$XG2 `basename $PWD`-st2-$i.pbs)
    else
	ST2=$ST2:$(qsub -W depend=afterany:$XG2 `basename $PWD`-st2-$i.pbs)
    fi
done

# Stage 3 - upper bound for veto algorithm
for i in `seq 1 $numScripts` ; do
   cp ../scripts/st3.pbs `basename $PWD`-st3-$i.pbs
   sed -i "s/nScripts=15/nScripts=$numScripts/g" `basename $PWD`-st3-$i.pbs
   sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-st3-$i.pbs
   qsub -W depend=afterany:$ST2 `basename $PWD`-st3-$i.pbs
done