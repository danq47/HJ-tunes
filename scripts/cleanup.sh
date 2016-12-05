basedir=$PBS_O_WORKDIR
cd $basedir

do_lhef=1
do_p8=1
do_NLOplots=1
do_libtunes=1

# First we want to merge all of the analysis files together

# Check that we have mergedata executables
if [ ! -e ../mergedata.f ] ; then
    echo;
    echo "Could not find mergedata.f file: Quitting."; echo
    exit
fi

rm -rf mergedata.f mergedata merge
cp ../mergedata.f .
gfortran -o mergedata mergedata.f

###########################
# First merge LHEF output #
###########################

if [ $do_lhef -eq 1 ] ; then
    cd LHEF-output
    cp ../mergedata .
    for i in `seq 1 7` ; do
	filesToCombine=`ls pwgLHEF_analysis-*-W$i.top`
	if [ "$filesToCombine" = "" ] ; then
	    echo "No LHEF W$i files found to merge."
	else
	    echo "About to merge files: "$filesToCombine

# Merge with both modes
	    for j in `seq 1 2` ; do
		echo "Output file: combined-output/pwgLHEF_analysis-W$i-combined-mode-$j.top"
		./mergedata $j $filesToCombine
		if [ -e fort.12 ] ; then
		    mv fort.12 ../combined-output/pwgLHEF_analysis-W$i-combined-mode-$j.top
		else
		    echo;
		    echo "Merge failed."
		    echo;
		fi
	    done
	fi
    done
# Finished merging, move back to base directory
    cd $basedir
fi

##############################
# Next, merge Pythia 8 files #
##############################

if [ $do_p8 -eq 1 ] ; then
    cd PYTHIA8-output
    cp ../mergedata .
    for i in `seq 1 7` ; do
        filesToCombine=`ls pwgPOWHEG+PYTHIA8-output-*-W$i.top`
        if [ "$filesToCombine" = "" ] ; then
            echo "No Pythia8 W$i files found to merge."
        else
            echo "About to merge files: "$filesToCombine

# Merge with both modes                                                                                  
            for j in `seq 1 2` ; do
                echo "Output file: combined-output/pwgPYTHIA8-W$i-combined-mode-$j.top"
                ./mergedata $j $filesToCombine
                if [ -e fort.12 ] ; then
                    mv fort.12 ../combined-output/pwgPYTHIA8-W$i-combined-mode-$j.top
                else
                    echo;
                    echo "Merge failed."
                    echo;
                fi
            done
        fi
    done
    cd $basedir
fi


#########################                                                                           
# Next, merge NLO files #                                                                           
#########################                                                                           

if [ $do_NLOplots -eq 1 ] ; then
    cd NLO-output
    cp ../mergedata .
    filesToCombine=`ls pwg-*-NLO.top`
    if [ "$filesToCombine" = "" ] ; then
        echo "No NLO files found to merge."
    else
        echo "About to merge files: "$filesToCombine
	
# Merge with both modes

        for j in `seq 1 2` ; do
            echo "Output file: combined-output/pwgNLO-combined-mode-$j.top"
            ./mergedata $j $filesToCombine
            if [ -e fort.12 ] ; then
                mv fort.12 ../combined-output/pwgNLO-combined-mode-$j.top
            else
                echo;
                echo "Merge failed."
                echo;
            fi
        done
    fi
    cd $basedir
fi

######################################################################
#                                                                    #
# Rescaling the output so that it is doubly differential in pT and Y #
#                                                                    #
######################################################################
if [ $do_libtunes -eq 1 ] ; then
    cd combined-output
    mkdir rescaled_output
    cd rescaled_output
    cp ../*.top .
    cp $basedir/rescale_top.py .
    for i in `ls *.top` ; do
	python rescale_top.py $i
    done
    rm *.top
    cd ..

####################################################
#                                                  #
# Outputting the rescaled results in libtunes form #
#                                                  #
####################################################

    cp -r rescaled_output libtunes_output
    cd libtunes_output
    cp $basedir/write_as_libtunes.py .
    for i in `ls *.top-rescaled` ; do
	python write_as_libtunes.py $i
    done
    rm *.top-rescaled
fi
cd $basedir

mkdir scripts
mv *pbs* scripts
mv *.sh.* scripts