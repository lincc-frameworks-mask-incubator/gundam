.. _outfiles:

============
Output Files
============

If ``write=True`` the code saves several files named as ``outfn`` plus an 
extension. According to the kind of correlation, these files are named as

+--------------------+---------------------------------------------------------+
| File extension     | Description                                             |
+====================+=========================================================+
| **.wrp,.xis,.wth** | Projected, redshift-space and angular correlation,      |
|                    | respectively (ASCII)                                    |
+--------------------+---------------------------------------------------------+
| **.dd,.rr,.dr**    | DD, RR, DR pair count arrays, in the corresponding      |
|                    | space bins (ASCII)                                      |
+--------------------+---------------------------------------------------------+
| **.cd,.cr**        | CD and CR pair count arrays, in the corresponding       |
|                    | space bins (ASCII)                                      |
+--------------------+---------------------------------------------------------+
| **.cnt**           | Pickled file with the **counts** dictionary object      |
|                    | (binary)                                                |
+--------------------+---------------------------------------------------------+
| **.par**           | Input parameters in JSON style dictionary (ASCII)       |
+--------------------+---------------------------------------------------------+

In order to properly capture log information from Python and Fortran code in
different environments, Gundam **always** generates two log files named as 
``outfn`` plus an extension as

+--------------------+---------------------------------------------------------+
| File extension     | Description                                             |
+====================+=========================================================+
| **.log**           | Log record produced by Python routines (ASCII)          |
+--------------------+---------------------------------------------------------+
| **.fortran.log**   | Log record produced by Fortran routines (ASCII)         |
+--------------------+---------------------------------------------------------+

If ``plot=True``, plots of correlation/counts are shown on screen and also 
saved to disk in eps format with extensions **.wrp.eps.**, **.xis.eps**, etc.
    
