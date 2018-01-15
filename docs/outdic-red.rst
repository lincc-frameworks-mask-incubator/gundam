==============
Redshift Space
==============

.. _outdicrcf:

Output Dictionary (rcf)
=======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dd, counts.rr, counts.xis, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| sl,sm,sr       | Left, mid, and right-points of redshift space bins (in Mpc) |
+----------------+-------------------------------------------------------------+
| xis,xiserr     | Redshift space correlation function and error               |
+----------------+-------------------------------------------------------------+
| dd             | DD pair count array in projected and radial bins            |
+----------------+-------------------------------------------------------------+
| rr             | RR pair count array in projected and radial bins            |
+----------------+-------------------------------------------------------------+
| dr             | DR pair count array in projected and radial bins            |
+----------------+-------------------------------------------------------------+
| bdd            | Boostrap DD pair count array in proj. and radial bins       |
+----------------+-------------------------------------------------------------+
| npt,npt1       | Nr. of points in data (D) and random sample (R), resp.      |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+


.. _outdicrccf:

Output Dictionary (rccf)
========================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.cd, counts.cr, counts.xis, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| sl,sm,sr       | Left, mid, and right-points of redshift space bins (in Mpc) |
+----------------+-------------------------------------------------------------+
| xis,xiserr     | Redshift space cross-correlation function and error         |
+----------------+-------------------------------------------------------------+
| cd             | CD pair count array in redshift space bins                  |
+----------------+-------------------------------------------------------------+
| cr             | CR pair count array in redshift space bins                  |
+----------------+-------------------------------------------------------------+
| bcd            | Boostrap CD pair count array in redshift space bins         |
+----------------+-------------------------------------------------------------+
| npt,npt1,npt2  | Nr. of points in samples D, R, and C, respectively          |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+


.. _outdicsA:

Output Dictionary (sA)
======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dd, counts.bd, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| sl,sm,sr       | Left, mid, and right-points of redshift space bins (in Mpc) |
+----------------+-------------------------------------------------------------+
| dd             | DD pair count array in redshift space bins                  |
+----------------+-------------------------------------------------------------+
| bdd            | Boostrap DD pair count array in redshift space bins         |
+----------------+-------------------------------------------------------------+
| npt            | Nr. of points in the sample (D)                             |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+


.. _outdicsC:

Output Dictionary (sC)
======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dr, counts.bdr, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| sl,sm,sr       | Left, mid, and right-points of redshift space bins (in Mpc) |
+----------------+-------------------------------------------------------------+
| dr             | DR pair count array in redshift space bins                  |
+----------------+-------------------------------------------------------------+
| bdr            | Boostrap DR pair count array in redshift space bins         |
+----------------+-------------------------------------------------------------+
| npt,npt1       | Nr. of points in the samples D and R, respectively          |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+
