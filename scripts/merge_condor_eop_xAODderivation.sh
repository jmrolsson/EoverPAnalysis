#!/bin/bash

cd $ROOTCOREBIN/..

echo "---> Merging data files:"
data_output=$(head -n 1 results/run_eop_data_xAODderivation_latest.log | sed -e 's/\/$//g')
echo hadd ${data_output}.root ${data_output}/fetch/hist-data15_13TeV_TrackCaloDxAOD-*.root
hadd ${data_output}.root ${data_output}/fetch/hist-data15_13TeV_TrackCaloDxAOD-*.root

echo "---> Merging MC files:"
mc_output=$(head -n 1 results/run_eop_mc_xAODderivation_latest.log | sed -e 's/\/$//g')
echo hadd ${mc_output}.root ${mc_output}/fetch/hist-mc15_13TeV_TrackCaloDxAOD-*.root
hadd ${mc_output}.root ${mc_output}/fetch/hist-mc15_13TeV_TrackCaloDxAOD-*.root

echo "---> Write to logfile:"
echo ${data_output}.root | sed -e 's/^results\///g' > results/merge_condor_eop_xAODderivation_latest.log
echo ${mc_output}.root | sed -e 's/^results\///g' >> results/merge_condor_eop_xAODderivation_latest.log

echo "--> All done!"
