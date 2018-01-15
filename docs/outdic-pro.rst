===============
Projected Space
===============

.. _outdicpcf:

Output Dictionary (pcf)
=======================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dd, counts.rr, counts.wrp, etc. It also 
stores the log and all input parameters.

+----------------+--------------------------------------------------------+
| Key(s)         | Description                                            |
+================+========================================================+
| rpl,rpm,rpr    | Left, mid, and right-points of projected bins (in Mpc) |
+----------------+--------------------------------------------------------+
| wrp,wrperr     | Projected correlation function and error               |
+----------------+--------------------------------------------------------+
| dd             | DD pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| rr             | RR pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| dr             | DR pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| bdd            | Boostrap DD pair count array in proj. and radial bins  |
+----------------+--------------------------------------------------------+
| npt,npt1       | Nr. of points in data (D) and random sample (R), resp. |
+----------------+--------------------------------------------------------+
| log            | Log record of Python routines                          |
+----------------+--------------------------------------------------------+
| logfortran     | Log record of Fortran routines                         |
+----------------+--------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)       |
+----------------+--------------------------------------------------------+


.. _outdicpccf:

Output Dictionary (pccf)
========================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.cd, counts.cr, counts.wrp, etc. It also 
stores the log and all input parameters.

+----------------+--------------------------------------------------------+
| Key(s)         | Description                                            |
+================+========================================================+
| rpl,rpm,rpr    | Left, mid, and right-points of projected bins (in Mpc) |
+----------------+--------------------------------------------------------+
| wrp,wrperr     | Projected cross-correlation function and error         |
+----------------+--------------------------------------------------------+
| cd             | CD pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| cr             | CR pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| bcd            | Boostrap CD pair count array in proj. and radial bins  |
+----------------+--------------------------------------------------------+
| npt,npt1,npt2  | Nr. of points in samples D, R, and C, respectively     |
+----------------+--------------------------------------------------------+
| log            | Log record of Python routines                          |
+----------------+--------------------------------------------------------+
| logfortran     | Log record of Fortran routines                         |
+----------------+--------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)       |
+----------------+--------------------------------------------------------+


.. _outdicrppiA:

Output Dictionary (rppiA)
=========================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dr, counts.bdr, counts.intpi, etc. It also 
stores the log and all input parameters.

+----------------+--------------------------------------------------------+
| Key(s)         | Description                                            |
+================+========================================================+
| rpl,rpm,rpr    | Left, mid, and right-points of projected bins (in Mpc) |
+----------------+--------------------------------------------------------+
| dd             | DD pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| bdd            | Boostrap DD pair count array in proj. and radial bins  |
+----------------+--------------------------------------------------------+
| intpi          | DD counts integrated along all radial bins             |
+----------------+--------------------------------------------------------+
| intpib         | Boostrap DD counts integrated along all radial bins    |
+----------------+--------------------------------------------------------+
| npt            | Nr. of points in the sample (D)                        |
+----------------+--------------------------------------------------------+
| log            | Log record of Python routines                          |
+----------------+--------------------------------------------------------+
| logfortran     | Log record of Fortran routines                         |
+----------------+--------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)       |
+----------------+--------------------------------------------------------+


.. _outdicrppiC:

Output Dictionary (rppiC)
=========================

Dictionary with attribute-style access, containing all counts and correlations 
accesible by field keys, e.g. counts.dr, counts.bdr, counts.intpi, etc. It also 
stores the log and all input parameters.

+----------------+--------------------------------------------------------+
| Key(s)         | Description                                            |
+================+========================================================+
| rpl,rpm,rpr    | Left, mid, and right-points of projected bins (in Mpc) |
+----------------+--------------------------------------------------------+
| dr             | DR pair count array in projected and radial bins       |
+----------------+--------------------------------------------------------+
| bdr            | Boostrap DR pair count array in proj. and radial bins  |
+----------------+--------------------------------------------------------+
| intpi          | DR counts integrated along all radial bins             |
+----------------+--------------------------------------------------------+
| intpib         | Boostrap DR counts integrated along all radial bins    |
+----------------+--------------------------------------------------------+
| npt,npt1       | Nr. of points in the samples D and R, respectively     |
+----------------+--------------------------------------------------------+
| log            | Log record of Python routines                          |
+----------------+--------------------------------------------------------+
| logfortran     | Log record of Fortran routines                         |
+----------------+--------------------------------------------------------+
| par            | Input parameters + runtime parameters (see xxxx)       |
+----------------+--------------------------------------------------------+
