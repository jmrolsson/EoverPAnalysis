Installing
==========

Checking out packages
---------------------

::

    setupATLAS
    lsetup 'rcsetup Base,X.Y.Z'
    git clone https://github.com/UCATLAS/xAODAnaHelpers
    git clone http://github.com/jmrolsson/EoverPAnalysis EoverPAnalysis

Compiling source
----------------

::

    rc clean
    rc find_packages
    rc compile
    rc make_par

.. note::
    - ``rc make_par`` creates a proof archive containing your RootCore packages, which is necessary if you want to run on a PROOF fram (e.g. condor).

More info
---------

  - For more information (about branches, svn, git, etc), please refer to the install page of |xAH|: `https://xaodanahelpers.readthedocs.io/en/master/Installing.html <https://xaodanahelpers.readthedocs.io/en/master/Installing.html>`_

  - RootCore twiki: `https://twiki.cern.ch/twiki/bin/view/AtlasComputing/RootCore <https://twiki.cern.ch/twiki/bin/view/AtlasComputing/RootCore>`_
