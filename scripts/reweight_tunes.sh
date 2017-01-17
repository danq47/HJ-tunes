# script to reweight all of our events with tunes reweighting and scale variation

if [ -e pwgevents.lhe ] ; then
    manyseeds=0
else
    manyseeds=1
fi

tarred_events=0 # 1 if events are in tarred format, 0 if not

if [ $tarred_events -eq 1 ] ; then # we will need to untar them
    for event_file in *.lhe.tar.gz ; do
	tar -xvzf $event_file 
    done
fi

###########################################################################
# Start a file to track how long things take.
> Timings-rwgt.txt
#
###########################################################################
              
if [ -e pwgevents-rwgt-*.lhe ] ; then                                                                                                                                        
    rm pwgevents-rwgt-*.lhe
fi

counter=0 # count which weight we are on

# First choose the POWHEG scales, we will have a separate weights group for each set of scales
for pwg_scales in {1..7}
do
    case $pwg_scales in
        1) facscfact=1.0 ; renscfact=1.0 ;;
        2) facscfact=0.5 ; renscfact=0.5 ;;
        3) facscfact=0.5 ; renscfact=1.0 ;;
        4) facscfact=1.0 ; renscfact=0.5 ;;
        5) facscfact=2.0 ; renscfact=1.0 ;;
        6) facscfact=1.0 ; renscfact=2.0 ;;
        7) facscfact=2.0 ; renscfact=2.0 ;;
    esac


# Turn on "tunes"
    cat  powheg.input-save | sed "s/parallelstage.*/parallelstage 4/ ; s/storeinfo_rwgt/compute_rwgt/ ; s/tunes.*/tunes 1/ " > powheg.input
    sed -i "s/facscfact .*/facscfact $facscfact/g" powheg.input
    sed -i "s/renscfact .*/renscfact $renscfact/g" powheg.input
    sed -i "s/lhrwgt_group_name .*/lhrwgt_group_name \'tunes reweight with POWHEG scales muR = $renscfact, muF = $facscfact\'/g" powheg.input


# Next, choose the MRT scales. Each set of POWHEG scale choices will have 7 MRT scale choices in the group
    for mrt_scales in {1..7}
    do
        case $mrt_scales in
            1) facscfact_mrt=1.0 ; renscfact_mrt=1.0 ;;
            2) facscfact_mrt=0.5 ; renscfact_mrt=0.5 ;;
            3) facscfact_mrt=0.5 ; renscfact_mrt=1.0 ;;
            4) facscfact_mrt=1.0 ; renscfact_mrt=0.5 ;;
            5) facscfact_mrt=2.0 ; renscfact_mrt=1.0 ;;
            6) facscfact_mrt=1.0 ; renscfact_mrt=2.0 ;;
            7) facscfact_mrt=2.0 ; renscfact_mrt=2.0 ;;
        esac
	
        sed -i "s/facscfact_mrt .*/facscfact_mrt $facscfact_mrt/g" powheg.input
        sed -i "s/renscfact_mrt .*/renscfact_mrt $renscfact_mrt/g" powheg.input
        sed -i "s/lhrwgt_id .*/lhrwgt_id \'muR_pwg = $renscfact, muF_pwg = $facscfact, muR_mrt = $renscfact_mrt, muF_mrt = $facscfact_mrt\'/g" powheg.input
        sed -i "s/lhrwgt_dscr .*/lhrwgt_dscr \'tunes reweighted with KRp=$muRp,KFp=$muFp,KRr=$muRr,KFr=$muFr\'/g" powheg.input
	
	counter=$((counter+1))


# Run the reweighting program with modified powheg.input files

    if [ $manyseeds -eq 0 ] ; then
        ./pwhg_main > run-rwgt.log
        mv pwgevents-rwgt.lhe pwgevents.lhe
    else

	    seed_list=`ls *.lhe | sed "s/pwgevents-//" | sed "s/.lhe//" `
	    for seed in $seed_list ; do
	        (echo -n ' ') >> Timings-rwgt.txt
	        (echo -n Reweighting event file $seed with weight number $counter of 49; date) >> Timings-rwgt.txt
	        ./pwhg_main <<EOF > run-rwgt-$seed.log 
$seed
pwgevents-$seed.lhe
EOF
	    
	       mv pwgevents-rwgt-$seed.lhe pwgevents-$seed.lhe

	    done

    fi

    done
done	

if [ -e log-files ] ; then
    mv run-rwgt-*.log log-files/
fi