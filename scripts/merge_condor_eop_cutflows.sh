#!/bin/bash

if [ $# -eq 0 ]

then
  echo "Usage: merge_condor_eop.sh run_condor_eop_logfile_name"

else

  cd $ROOTCOREBIN/..

  logfile=$1
  echo
  echo "Reading from logfile: "${logfile}
  echo

  count=0
  while read line; do
    if [ $((count%2)) -eq 0 ]
    then
      file_tag=$(echo ${line} | sed -e 's/^.*\///g' | sed -e 's/\.txt$//g')
      echo "---> Processing output of: "${file_tag}
    else
      condor_output=${line}
      echo "condor output files: "${condor_output}
      echo "merging files:"
      echo hadd -f ${condor_output}-cutflows.root ${condor_output}/fetch/data-cutflow/hist-${file_tag}-*.root
      echo "writing to logfile:"
      if [ $((count)) -eq 1 ]
      then
        echo ${condor_output}.root | sed -e 's/^results\///g' > results/merge_condor_eop_lowmu_latest.log
      else
        echo ${condor_output}.root | sed -e 's/^results\///g' >> results/merge_condor_eop_lowmu_latest.log
      fi
    fi
    ((count++));
  done <${logfile};

  echo "--> All done!"

fi
