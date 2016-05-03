#!/bin/bash
if [ $# -eq 0 ]

  then
    echo "Usage: source run_eop_xAODderivation_condor.sh tag"

  else

    cd $ROOTCOREBIN/..

    tag=$1
    today=$(date +"%Y%m%d")

    # echo "---> Running data:"
    # echo xAH_run.py --files EoverP/filelists/data15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_xAODderivation_data.py --submitDir results/eop_data_xAODderivation_condor_${today}_${tag} --verbose --force condor --optFilesPerWorker 80
    # xAH_run.py --files EoverP/filelists/data15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_xAODderivation_data.py --submitDir results/eop_data_xAODderivation_condor_${today}_${tag} --verbose --force condor --optFilesPerWorker 80
    #
    # echo "---> Running MC:"
    # echo xAH_run.py --files EoverP/filelists/mc15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_xAODderivation_mc.py --submitDir results/eop_mc_xAODderivation_condor_${today}_${tag} --verbose --force condor --optFilesPerWorker 20
    # xAH_run.py --files EoverP/filelists/mc15_13TeV_TrackCaloDxAOD.txt --inputList --config EoverP/scripts/config_eop_xAODderivation_mc.py --submitDir results/eop_mc_xAODderivation_condor_${today}_${tag} --verbose --force condor --optFilesPerWorker 20
    #
    # echo "---> Merging data files:"
    # cd results/eop_data_xAODderivation_condor_${today}_${tag}/fetch/
    # echo hadd ../../data15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root hist-data15_13TeV_TrackCaloDxAOD-*.root > /dev/null
    # hadd ../../data15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root hist-data15_13TeV_TrackCaloDxAOD-*.root > /dev/null
    # cd -

    echo "---> Merging MC files:"
    cd results/eop_mc_xAODderivation_condor_${today}_${tag}/fetch/
    echo hadd ../../mc15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root hist-mc15_13TeV_TrackCaloDxAOD-*.root > /dev/null
    hadd ../../mc15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root hist-mc15_13TeV_TrackCaloDxAOD-*.root > /dev/null
    cd -

    echo "---> Write to logfile:"
    echo  data15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root > results/run_eop_xAODderivation_condor_${today}_${tag}.log

    echo  mc15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root >> results/run_eop_xAODderivation_condor_${today}_${tag}.log
    echo  data15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root > results/run_eop_xAODderivation_condor_latest.log

    echo  mc15_13TeV_TrackCaloDxAOD_condor_${today}_${tag}.root >> results/run_eop_xAODderivation_condor_latest.log

    echo "--> All done!"

fi
