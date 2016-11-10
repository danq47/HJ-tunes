# Stage 3 - upper bound for veto algorithm                                                                                         

numScripts=25
 
for i in `seq 1 $numScripts` ; do                                                                                                  
cp ../scripts/st3.pbs `basename $PWD`-st3-$i.pbs                                                                               
sed -i "s/nScripts=15/nScripts=$numScripts/g" `basename $PWD`-st3-$i.pbs                                                       
sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-st3-$i.pbs                                                        
qsub `basename $PWD`-st3-$i.pbs                                                                      
done