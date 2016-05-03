#!/bin/bash
if [ $# -eq 0 ]

  then
    echo "Usage: source run_eop_xAODderivation_test.sh tag"

  else

    tag=$1
    today=$(date +"%Y%m%d")

    echo "Running data:"
    echo xAH_run.py --files EoverP/filelists/data15_13TeV_TrackCaloDxAOD_test1.txt --inputList --config EoverP/scripts/config_eop_xAODderivation.py --submitDir results/eop_data_xAODderivation_test_${today}_${tag} --verbose --force direct
    xAH_run.py --files EoverP/filelists/data15_13TeV_TrackCaloDxAOD_test1.txt --inputList --config EoverP/scripts/config_eop_xAODderivation.py --submitDir results/eop_data_xAODderivation_test_${today}_${tag} --verbose --force direct

    echo "Running MC:"
    echo xAH_run.py --files EoverP/filelists/mc15_13TeV_TrackCaloDxAOD_test1.txt --inputList --config EoverP/scripts/config_eop_xAODderivation.py --submitDir results/eop_mc_xAODderivation_test_${today}_${tag} --verbose --force direct
    xAH_run.py --files EoverP/filelists/mc15_13TeV_TrackCaloDxAOD_test1.txt --inputList --config EoverP/scripts/config_eop_xAODderivation.py --submitDir results/eop_mc_xAODderivation_test_${today}_${tag} --verbose --force direct

    # Write to logfile
    echo eop_data_xAODderivation_test_${today}_${tag} > results/run_eop_xAODderivation_test_${today}_${tag}.log
    echo eop_mc_xAODderivation_test_${today}_${tag} >> results/run_eop_xAODderivation_test_${today}_${tag}.log
    echo eop_data_xAODderivation_test_${today}_${tag} > results/run_eop_xAODderivation_test_latest.log
    echo eop_mc_xAODderivation_test_${today}_${tag} >> results/run_eop_xAODderivation_test_latest.log

fi
