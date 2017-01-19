#!/bin/bash
if [ $# -eq 0 ]

  then
    echo "Usage: source run_condor_all_eop_lowmu_runII_general.sh tag"

  else

    cd $ROOTCOREBIN/..

    tag=$1
    today=$(date +"%Y%m%d")
    files_data=EoverPAnalysis/filelists/data15_13TeV_lowmu_all_NEW.txt
    files_mc=EoverPAnalysis/filelists/mc15_13TeV_lowmu_all_NEW.txt

    mkdir -p results

    echo "---> Running data:"
    echo xAH_run.py --files ${files_data} --inputList --config EoverPAnalysis/scripts/config_eop_data_lowmu_runII_general.py --submitDir results/condor_all_eop_lowmu_runII_general_data_${today}_${tag} --verbose --force condor --optFilesPerWorker 40
    xAH_run.py --files ${files_data} --inputList --config EoverPAnalysis/scripts/config_eop_data_lowmu_runII_general.py --submitDir results/condor_all_eop_lowmu_runII_general_data_${today}_${tag} --verbose --force condor --optFilesPerWorker 40

    echo "---> Running MC:"
    echo xAH_run.py --files ${files_mc} --inputList --config EoverPAnalysis/scripts/config_eop_mc_lowmu_runII_general.py --submitDir results/condor_all_eop_lowmu_runII_general_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 20
    xAH_run.py --files ${files_mc} --inputList --config EoverPAnalysis/scripts/config_eop_mc_lowmu_runII_general.py --submitDir results/condor_all_eop_lowmu_runII_general_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 20

    echo "---> Write to logfile:"
    echo ${files_data} > results/run_condor_eop_lowmu_runII_general.log
    echo results/condor_all_eop_lowmu_runII_general_data_${today}_${tag} >> results/run_condor_eop_lowmu_runII_general.log
    echo ${files_mc} >> results/run_condor_eop_lowmu_runII_general.log
    echo results/condor_all_eop_lowmu_runII_general_mc_${today}_${tag} >> results/run_condor_eop_lowmu_runII_general.log

    echo "--> Jobs submitted!"
    echo "source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.sh $ROOTCOREBIN/../results/run_condor_eop_lowmu_runII_general.log # when condor jobs are finished to merge output files"

fi
