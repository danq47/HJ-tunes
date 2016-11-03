startingSeed=1

cp ../scripts/DQ-0.pbs `basename $PWD`-initialisation.pbs
sed -i "s/startSeed=1/startSeed=$startingSeed/g" `basename $PWD`-initialisation.pbs

# The jobs that follow need the initialization stage to precede them
FIRST=$(qsub `basename $PWD`-initialisation.pbs)

for i in `seq 1 15` ; do

cp ../scripts/DQ-$i.pbs `basename $PWD`-$i.pbs
sed -i "s/scriptNumber-1)\ )\ +\ 1/scriptNumber-1)\ )\ +\ $startingSeed/g" `basename $PWD`-$i.pbs
qsub -W depend=afterany:$FIRST `basename $PWD`-$i.pbs

done