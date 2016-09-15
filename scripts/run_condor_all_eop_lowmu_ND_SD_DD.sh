#!/bin/bash
if [ $# -eq 0 ]

  then
    echo "Usage: source run_condor_all_eop_lowmu_ND_SD_DD.sh tag"

  else

    cd $ROOTCOREBIN/..

    tag=$1
    today=$(date +"%Y%m%d")
    files_ND=EoverPAnalysis/filelists/mc15_13TeV_lowmu_ND.txt
    files_SD=EoverPAnalysis/filelists/mc15_13TeV_lowmu_SD.txt
    files_DD=EoverPAnalysis/filelists/mc15_13TeV_lowmu_DD.txt

    mkdir -p results

    echo "---> Running ND, SD, and DD lowmu samples:"
    echo xAH_run.py --files ${files_ND} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_ND_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 30
    xAH_run.py --files ${files_ND} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_ND_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 30
    echo xAH_run.py --files ${files_SD} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_SD_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 30
    xAH_run.py --files ${files_SD} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_SD_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 30
    echo xAH_run.py --files ${files_DD} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_DD_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 30
    xAH_run.py --files ${files_DD} --inputList --config EoverPAnalysis/scripts/config_eop_mc.py --submitDir results/condor_all_eop_DD_mc_${today}_${tag} --verbose --force condor --optFilesPerWorker 30

    echo "---> Write to logfile:"
    # echo "# ---> "$(date +"%Y-%m-%d:%H:%M:%S") >> results/run_condor_eop_lowmu_ND_SD_DD.log
    # echo ${files_ND} >> results/run_condor_eop_lowmu.log
    # echo results/condor_all_eop_ND_mc_${today}_${tag} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    # echo ${files_SD} >> results/run_condor_eop_lowmu.log
    # echo results/condor_all_eop_SD_mc_${today}_${tag} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    # echo ${files_DD} >> results/run_condor_eop_lowmu.log
    # echo results/condor_all_eop_DD_mc_${today}_${tag} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    echo ${files_ND} > results/run_condor_eop_lowmu_ND_SD_DD.log
    echo results/condor_all_eop_ND_mc_${today}_${tag} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    echo ${files_SD} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    echo results/condor_all_eop_SD_mc_${today}_${tag} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    echo ${files_DD} >> results/run_condor_eop_lowmu_ND_SD_DD.log
    echo results/condor_all_eop_DD_mc_${today}_${tag} >> results/run_condor_eop_lowmu_ND_SD_DD.log

    echo "--> Jobs submitted!"
    echo "source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.sh $ROOTCOREBIN/../results/run_condor_eop_lowmu_ND_SD_DD.log # when condor jobs are finished to merge output files"

fi
