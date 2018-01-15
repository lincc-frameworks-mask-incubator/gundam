=============
Angular Space
=============

.. _outdicacf:

Output Dictionary (acf)
=======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dd, counts.rr, counts.wth, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| thl,thm,thr    | Left, mid, and right-points of anglar bins (in deg)         |
+----------------+-------------------------------------------------------------+
| wth,wtherr     | Angular correlation function and error                      |
+----------------+-------------------------------------------------------------+
| dd             | DD pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| rr             | RR pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| dr             | DR pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| bdd            | Boostrap DD pair count array in angular bins                |
+----------------+-------------------------------------------------------------+
| npt,npt1       | Nr. of points in data (D) and random sample (R), resp.      |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+


.. _outdicaccf:

Output Dictionary (accf)
========================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.cd, counts.cr, counts.wth, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| thl,thm,thr    | Left, mid, and right-points of angular bins (in Mpc)        |
+----------------+-------------------------------------------------------------+
| wht,wtherr     | Angular cross-correlation function and error                |
+----------------+-------------------------------------------------------------+
| cd             | CD pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| cr             | CR pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| bcd            | Boostrap CD pair count array in angular bins                |
+----------------+-------------------------------------------------------------+
| npt,npt1,npt2  | Nr. of points in samples D, R, and C, respectively          |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+


.. _outdicthA:

Output Dictionary (thA)
=======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dd, counts.bdd, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| thl,thm,thr    | Left, mid, and right-points of angular bins (deg)           |
+----------------+-------------------------------------------------------------+
| dd             | DD pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| bdd            | Boostrap DD pair count array in angular bins                |
+----------------+-------------------------------------------------------------+
| npt            | Nr. of points in the sample (D)                             |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+


.. _outdicthC:

Output Dictionary (thC)
=======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dr, counts.bdr, etc. It also 
stores the log and all input parameters.

+----------------+-------------------------------------------------------------+
| Key(s)         | Description                                                 |
+================+=============================================================+
| thl,thm,thr    | Left, mid, and right-points of angular bins (deg)           |
+----------------+-------------------------------------------------------------+
| dr             | DR pair count array in angular bins                         |
+----------------+-------------------------------------------------------------+
| bdr            | Boostrap DR pair count array in angular bins                |
+----------------+-------------------------------------------------------------+
| npt,npt1       | Nr. of points in the samples D and R, respectively          |
+----------------+-------------------------------------------------------------+
| log            | Log record of Python routines                               |
+----------------+-------------------------------------------------------------+
| logfortran     | Log record of Fortran routines                              |
+----------------+-------------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)            |
+----------------+-------------------------------------------------------------+
