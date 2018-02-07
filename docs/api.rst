API Reference
=============

Main Routines (projected space)
-------------------------------

**These are the main routines for counting pairs and estimating correlation functions.
They encapsulate all required steps to provide a true end-user experience : you
only need to provide data+parameters to receive final results, and perhaps that
amazing plot of your next paper**

.. autofunction:: gundam.pcf
.. autofunction:: gundam.pccf
.. autofunction:: gundam.rppi_A
.. autofunction:: gundam.rppi_C


Main Routines (redshift space)
------------------------------

**These are the main routines for counting pairs and estimating correlation functions.
They encapsulate all required steps to provide a true end-user experience : you
only need to provide data+parameters to receive final results, and even that
amazing plot of your next paper**

.. autofunction:: gundam.rcf
.. autofunction:: gundam.rccf
.. autofunction:: gundam.s_A
.. autofunction:: gundam.s_C


Main Routines (angular space)
-----------------------------

**These are the main routines for counting pairs and estimating correlation functions.
They encapsulate all required steps to provide a true end-user experience : you
only need to provide data+parameters to receive final results, and even that
amazing plot of your next paper**

.. autofunction:: gundam.acf
.. autofunction:: gundam.accf
.. autofunction:: gundam.th_A
.. autofunction:: gundam.th_C


Auxiliary Routines
------------------

**These are helper routines that make it easy to work with correlations, visualize
count data, optimize grid size, sort data, etc. For examples of how to use them
check their individual help, or even better, the source code of** :func:`gundam.pcf`

.. autofunction:: gundam.bestSKgrid2d
.. autofunction:: gundam.bestSKgrid3d
.. autofunction:: gundam.cfindex
.. autofunction:: gundam.cnttable
.. autofunction:: gundam.makebins
.. autofunction:: gundam.pixsort
.. autofunction:: gundam.qprint
.. autofunction:: gundam.radec2xyz
.. autofunction:: gundam.tpcf
.. autofunction:: gundam.tpccf
.. autofunction:: gundam.tpcf_wrp
.. autofunction:: gundam.tpccf_wrp


Plotting Routines
-----------------

.. autofunction:: gundam.cntplot
.. autofunction:: gundam.cntplot2D
.. autofunction:: gundam.comparecf
.. autofunction:: gundam.fitpowerlaw
.. autofunction:: gundam.plotcf


Comprehensive Auxiliary Routines
--------------------------------
**These are helper routines that perform many common tasks (in the Python side)
during a correlation run, such as collecting counts, i/o functions, initialization,
logging, etc. For a nice example of how to use them, see source code of**
:func:`gundam.pcf`

.. autofunction:: gundam.addlog
.. autofunction:: gundam.aggcounts
.. autofunction:: gundam.aggcountsb
.. autofunction:: gundam.allequal
.. autofunction:: gundam.addpvcols
.. autofunction:: gundam.bound2d
.. autofunction:: gundam.bound3d
.. autofunction:: gundam.cross0guess
.. autofunction:: gundam.buildoutput
.. autofunction:: gundam.buildoutputC
.. autoclass:: gundam.capture
.. autofunction:: gundam.check_kind
.. autofunction:: gundam.closelog
.. autofunction:: gundam.finalize
.. autofunction:: gundam.init
.. autofunction:: gundam.initialize
.. autofunction:: gundam.logcallinfo
.. autofunction:: gundam.logtimming
.. autofunction:: gundam.readcounts
.. autofunction:: gundam.readpars
.. autofunction:: gundam.savecounts
.. autofunction:: gundam.savepars
.. autofunction:: gundam.setlogs
.. autofunction:: gundam.writeasc_cf
.. autofunction:: gundam.writeasc_rppicounts
.. autofunction:: gundam.writeasc_counts
.. autofunction:: gundam.watch

.. _api_fwrappers:

Fortran Wrapper Routines
------------------------

These are the wrappers that call Fortran pair counting routines. They are intended
to: (1) choose the fastest counting routine depending whether weights, 
bootstrap errors, etc. are requested or not; and (2): concentrate all mayor 
Fortran calling in one place

Unless you are building your own custom script, heavily modifying, or extending
the core algorithms, you normally should not need to use these functions

.. autofunction:: gundam.pairs_auto
.. autofunction:: gundam.pairs_cross
