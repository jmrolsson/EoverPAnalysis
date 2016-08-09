# EoverPAnalysis

<a href="http://eoverp.readthedocs.io/en/latest/">The documentation is hosted on ReadTheDocs.org</a> 

<h2>Basic usage:</h2>

<h3>Setup</h3>
<p>mkdir myAnalysis; cd myAnalysis</p>
<p>git clone http://github.com/UCATLAS/xAODAnaHelpers xAODAnaHelpers</p>
<p>git clone http://github.com/jmrolsson/EoverPAnalysis EoverPAnalysis</p>
<p>lsetup 'rcsetup Base,2.4.X'</p>
<p>rc clean && rc find_packages && rc compile && rc make_par</p>

<h3>Test run</h3>
<p>mkdir results</p>
<p>xAH_run.py --files EoverPAnalysis/filelists/data15_13TeV_lowmu_test1.txt --inputList --config EoverPAnalysis/scripts/config_eop_data.py --submitDir results/eop_data_test_0 --verbose --force direct</p>
