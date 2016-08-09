# EoverPAnalysis
<a href="http://eoverp.readthedocs.io/en/latest/">The complete documentation is hosted on ReadTheDocs.org</a> 

For questions please contact: joakim.olsson[at]cern.ch

## Setup

```
mkdir myAnalysis; cd myAnalysis
git clone http://github.com/UCATLAS/xAODAnaHelpers xAODAnaHelpers
git clone http://github.com/jmrolsson/EoverPAnalysis EoverPAnalysis
lsetup 'rcsetup Base,2.4.X'
rc clean && rc find_packages && rc compile && rc make_par
```

## Running

### Local test run

```
mkdir results
xAH_run.py --files EoverPAnalysis/filelists/data15_13TeV_lowmu_test1.txt --inputList --config EoverPAnalysis/scripts/config_eop_data_lowmu.py --submitDir results/eop_data_test_0 --verbose --force direct
```

### First condor test run

```
source EoverPAnalysis/scripts/run_condor_test_eop_lowmu.sh 0 # where '0' is a tag for the run
```

The output will then be located in 'results', e.g. results/condor_test_eop_lowmu_{mc,data}_YYYYMMDD_0/</p>

### The condor output histograms and cutflows can easily be merged, just run (after condor has finished)

```
source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.py $ROOTCOREBIN/../results/run_condor_eop_lowmu_latest.log
```

## Configuration

### config_* scripts

In 'scripts' you'll find files with names like 'config_*' (ex. 'config_data.py'). These files set the run options, i.e. what event and track selection to apply, what histograms to make, etc. Create your own as needed! 

### run_condor_* scripts

In 'scripts' you'll also find files with names such as 'run_condor_*' (ex. 'run_condor_test_eop_lowmu.sh'). These let you automate submission to condor (see 'First condor test run' above for instructions).
