#!/bin/bash
if [ $# -eq 0 ]

  then
    echo "Usage: source run_condor_all_eop_xAODderivation.sh tag"

  else

    cd $ROOTCOREBIN/..

    tag=$1
    today=$(date +"%Y%m%d")

    echo "---> Running data:"
    echo xAH_run.py --files EoverP/filelists/data15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_data_xAODderivation.py --submitDir results/condor_all_eop_data_xAODderivation_${today}_${tag} --verbose --force condor --optFilesPerWorker 20
    xAH_run.py --files EoverP/filelists/data15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_data_xAODderivation.py --submitDir results/condor_all_eop_data_xAODderivation_${today}_${tag} --verbose --force condor --optFilesPerWorker 20

    echo "---> Running MC:"
    echo xAH_run.py --files EoverP/filelists/mc15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_mc_xAODderivation.py --submitDir results/condor_all_eop_mc_xAODderivation_${today}_${tag} --verbose --force condor --optFilesPerWorker 20
    xAH_run.py --files EoverP/filelists/mc15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_mc_xAODderivation.py --submitDir results/condor_all_eop_mc_xAODderivation_${today}_${tag} --verbose --force condor --optFilesPerWorker 20

    echo "---> Write to logfile:"
    # logfiles of all runs
    echo results/condor_all_eop_data_xAODderivation_${today}_${tag} >> results/run_condor_all_eop_xAODderivation.log
    echo results/condor_all_eop_mc_xAODderivation_${today}_${tag} >> results/run_condor_all_eop_xAODderivation.log
    # logfile read by ./scripts/merge_condor_eop_xAODderivation.py
    echo results/condor_all_eop_data_xAODderivation_${today}_${tag} > results/run_eop_data_xAODderivation_latest.log
    echo results/condor_all_eop_mc_xAODderivation_${today}_${tag} > results/run_eop_mc_xAODderivation_latest.log

    echo "--> Jobs submitted!"
    echo "run ./scripts/merge_condor_eop_xAODderivation.py when  condor jobs are finished to merge output"

fi
