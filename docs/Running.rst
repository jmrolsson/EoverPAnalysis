Running
========

Grid proxy
----------

If your datasets are located on the grid (the ones in the file lists that comes with this package are), you need to have a valid grid proxy in order to access them.

::

    voms-proxy-init -voms atlas

If you haven't done so already, you might want to add the following lines to your ~/.bash_profile:

::

    alias grid="voms-proxy-init -voms atlas -out $HOME/.globus/gridproxy.cert -valid 1000:00"
    export X509_USER_PROXY=$HOME/.globus/gridproxy.cert

The datasets that come with the default package are located on MWT2_UC_LOCALGROUPDISK, and are accessed via FAX. To set up fax, do:

::

    lsetup fax; fax-get-best-redirector

Local test run
--------------

::

    mkdir results
    xAH_run.py --files $ROOTCOREBIN/../EoverPAnalysis/filelists/data15_13TeV_lowmu_test1.txt --inputList --config $ROOTCOREBIN/../EoverPAnalysis/scripts/config_eop_data_lowmu.py --submitDir $ROOTCOREBIN/../results/eop_data_test_0 --verbose --force direct

First condor test run
---------------------

::

    source $ROOTCOREBIN/../EoverPAnalysis/scripts/run_condor_test_eop_lowmu.sh 0 # where '0' is a tag for the run

.. note::
    The output will then be located in 'results', e.g. $ROOTCOREBIN/../results/condor_test_eop_lowmu_{mc,data}_YYYYMMDD_0/

The condor output histograms and cutflows can easily be merged, just run (after condor has finished)

::
    
    source $ROOTCOREBIN/../EoverPAnalysis/scripts/merge_condor_eop.py $ROOTCOREBIN/../results/run_condor_eop_lowmu_latest.log
