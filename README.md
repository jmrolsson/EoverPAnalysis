# EoverPAnalysis

<a href="http://eoverp.readthedocs.io/en/latest/">The complete documentation is hosted on ReadTheDocs.org</a> 

<h2>Setup</h2>
<p>mkdir myAnalysis; cd myAnalysis</p>
<p>git clone http://github.com/UCATLAS/xAODAnaHelpers xAODAnaHelpers</p>
<p>git clone http://github.com/jmrolsson/EoverPAnalysis EoverPAnalysis</p>
<p>lsetup 'rcsetup Base,2.4.X'</p>
<p>rc clean && rc find_packages && rc compile && rc make_par</p>

<h2>Running:</h2>

<h3>Local test run</h3>
<p>mkdir results</p>
<p>xAH_run.py --files EoverPAnalysis/filelists/data15_13TeV_lowmu_test1.txt --inputList --config EoverPAnalysis/scripts/config_eop_data_lowmu.py --submitDir results/eop_data_test_0 --verbose --force direct</p>

<h3>First condor test run</h3>
<p>source EoverPAnalysis/scripts/run_condor_test_eop_lowmu.sh 0 # where '0' is a tag for the run</p>
<h4>The output will then be located in 'results':</h4>
<p>E.g. results/condor_test_eop_lowmu_{mc,data}_YYYYMMDD_0/</p>
<h4>The condor output histograms and cutflows can easily be merged, just run (after condor has finished):</h4> 
<p>source $RootCoreBin/../EoverPAnalysis/scripts/merge_condor_eop.py $RootCoreBin/../results/run_condor_eop_lowmu_latest.log

<h2>Configuration</h2>

<h3>config_* scripts</h3>

<p>In 'scripts' you'll find files with names like 'config_*' (ex. 'config_data.py'). These files set the run options, i.e. what event and track selection to apply, what histograms to make, etc. Create your own as needed!</p>

<h3>run_condor_* scripts</h3>

<p>In 'scripts' you'll also find files with names such as 'run_condor_*' (ex. 'run_condor_test_eop_lowmu.sh'). These let you automate submission to condor (see 'First condor test run' above for instructions).</p>

For questions please contact: joakim.olsson[AT]cern.ch
