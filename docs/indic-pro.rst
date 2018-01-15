===============
Projected Space
===============

.. _indic-pcf:

Input Parameters Dictionary (pcf)
=================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('pcf')                              |
+-------------+-------------------------------------------------------------------+
| h0          | Hubble constant [km/s/Mpc]. Default=100.                          |
+-------------+-------------------------------------------------------------------+
| omegam      | Omega matter. Default=0.3                                         |
+-------------+-------------------------------------------------------------------+
| omegal      | Omega lambda. Default=0.7                                         |
+-------------+-------------------------------------------------------------------+
| autogrid    | If ``True`` choose the optimum nr. of cells                       |
|             | ``(mxh1,mxh2,mxh3)`` of the skip table (SK). Default=True         |
+-------------+-------------------------------------------------------------------+
| dens        | Custom nr. of particles per SK cell used when ``autogrid=True``.  |
|             | No need to specify unless desired. Default=None                   |
+-------------+-------------------------------------------------------------------+
| mxh1        | Nr. of DEC cells of the SK table. Only relevant if                |
|             | ``autogrid=False``. Default=30                                    |
+-------------+-------------------------------------------------------------------+
| mxh2        | Nr. of RA cells of the SK table. Only relevant if                 |
|             | ``autogrid=False``. Default=180                                   |
+-------------+-------------------------------------------------------------------+
| mxh3        | Nr. of DCOM cells of the SK table. Only relevant if               |
|             | ``autogrid=False``. See :ref:`notemxh3` below. Default=40         |
+-------------+-------------------------------------------------------------------+
| pxorder     | Pixel ordering method. See :ref:`pxorder` for details and         |
|             | options. Default='natural'                                        |
+-------------+-------------------------------------------------------------------+
| doboot      | If ``True``, calculate bootstrap counts and error bars.           |
|             | Default=False                                                     |
+-------------+-------------------------------------------------------------------+
| nbts        | Nr. of bootstrap samples. Only relevant if ``doboot=True``.       |
|             | Default=50                                                        |
+-------------+-------------------------------------------------------------------+
| bseed       | Fixed seed for boostrap RNG. Always set ``bseed>0`` if running    |
|             | in paralell. Default=12345                                        |
+-------------+-------------------------------------------------------------------+
| wfib        | If ``True`` apply SDSS fiber correction for pairs closer than     |
|             | 55". See *wfiber()* Fortran function. Default=False               |
+-------------+-------------------------------------------------------------------+
| nsepp       | Nr. of projected separation bins. Default=22                      |
+-------------+-------------------------------------------------------------------+
| seppmin     | Minimum projected distance to consider [Mpc/h]. Default=0.01      |
+-------------+-------------------------------------------------------------------+
| dsepp       | Size of projected bins (in dex if log bins). Default=0.15         |
+-------------+-------------------------------------------------------------------+
| logsepp     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
+-------------+-------------------------------------------------------------------+
| nsepv       | Nr. of radial separation bins. Default=1                          |
+-------------+-------------------------------------------------------------------+
| dsepv       | Size of radial bins [Mpc/h]. Default=40.                          |
+-------------+-------------------------------------------------------------------+
| calcdist    | If ``False``, take comov. distances from input tables             |
|             | instead of calculating them. Default=True                         |
+-------------+-------------------------------------------------------------------+
| file        | File name of **data** sample. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| file1       | File name of **random** sample. Only informative. Default=''      |
+-------------+-------------------------------------------------------------------+
| description | Short description of the run. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| estimator   | Statistical estimator of the correlation function.                |
|             | Default='NAT' |br|                                                |
|             | * 'NAT' : Natural estimator -> :math:`DD/RR-1` |br|               |
|             | * 'HAM' : Hamilton estimator -> :math:`DD*RR/DR^{2}-1` |br|       |
|             | * 'LS' : Landy-Szalay estimator -> :math:`(DD-2DR+RR)/RR` |br|    |
|             | * 'DP' : Davis-Peebles estimator -> :math:`DD/DR-1` |br|          |
+-------------+-------------------------------------------------------------------+
| cra,cdec,   | Column names in **data** sample table (tab).                      |
| cred,cwei,  | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom       |                                                                   |
+-------------+-------------------------------------------------------------------+
| cra1,cdec1, | Column names in **random** sample (tab1).                         |
| cred1,cwei1 | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom1      |                                                                   |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _indic-pccf:

Input Parameters Dictionary (pccf)
==================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('pccf')                             |
+-------------+-------------------------------------------------------------------+
| h0          | Hubble constant [km/s/Mpc]. Default=100.                          |
+-------------+-------------------------------------------------------------------+
| omegam      | Omega matter. Default=0.3                                         |
+-------------+-------------------------------------------------------------------+
| omegal      | Omega lambda. Default=0.7                                         |
+-------------+-------------------------------------------------------------------+
| autogrid    | If ``autogrid=True`` choose the optimum nr. of cells              |            
|             | ``(mxh1,mxh2,mxh3)`` of the skip table (SK). Default=True         |
+-------------+-------------------------------------------------------------------+
| dens        | Custom nr. of particles per SK cell used when ``autogrid=True``.  |
|             | No need to specify unless desired. Default=None                   |
+-------------+-------------------------------------------------------------------+
| mxh1        | Nr. of DEC cells of the SK table. Only relevant if                |
|             | ``autogrid=False``. Default=30                                    |
+-------------+-------------------------------------------------------------------+
| mxh2        | Nr. of RA cells of the SK table. Only relevant if                 |
|             | ``autogrid=False``. Default=180                                   |
+-------------+-------------------------------------------------------------------+
| mxh3        | Nr. of DCOM cells of the SK table. Only relevant if               |
|             | ``autogrid=False``. See :ref:`notemxh3` below. Default=40         |
+-------------+-------------------------------------------------------------------+
| doboot      | If ``True``, calculate bootstrap counts and error bars.           |
|             | Default=False                                                     |
+-------------+-------------------------------------------------------------------+
| nbts        | Nr. of bootstrap samples. Only relevant if ``doboot=True``.       |
|             | Default=50                                                        |
+-------------+-------------------------------------------------------------------+
| bseed       | Fixed seed for boostrap RNG. Always set ``bseed>0`` if running    |
|             | in paralell. Default=1245                                         |
+-------------+-------------------------------------------------------------------+
| wfib        | If ``True`` apply SDSS fiber correction for pairs closer than     |
|             | 55". See *wfiber()* Fortran function. Default=False               |
+-------------+-------------------------------------------------------------------+
| nsepp       | Nr. of projected separation bins. Default=22                      |
+-------------+-------------------------------------------------------------------+
| seppmin     | Minimum projected distance to consider [Mpc/h]. Default=0.01      |
+-------------+-------------------------------------------------------------------+
| dsepp       | Size of projected bins (in dex if log bins). Default=0.15         |
+-------------+-------------------------------------------------------------------+
| logsepp     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
+-------------+-------------------------------------------------------------------+
| nsepv       | Nr. of radial separation bins. Default=1                          |
+-------------+-------------------------------------------------------------------+
| dsepv       | Size of radial bins [Mpc/h]. Default=40.                          |
+-------------+-------------------------------------------------------------------+
| calcdist    | If ``False``, take comov. distances from input tables             |
|             | instead of calculating them. Default=True                         |
+-------------+-------------------------------------------------------------------+
| file        | File name of **data** sample. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| file1       | File name of **random** sample. Only informative. Default=''      |
+-------------+-------------------------------------------------------------------+
| file2       | File name of **cross** sample. Only informative. Default=''       |
+-------------+-------------------------------------------------------------------+
| description | Short description of the run. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| estimator   | Statistical estimator of the correlation function.                |
|             | Default='DP' |br|                                                 |
|             | * 'DP' : Davis-Peebles estimator -> :math:`CD/CR - 1`             |
+-------------+-------------------------------------------------------------------+
| cra,cdec,   | Column names in **data** sample table (tab).                      |
| cred,cwei,  | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom       |                                                                   |
+-------------+-------------------------------------------------------------------+
| cra1,cdec1, | Column names in **random** sample (tab1).                         |
| cred1,cwei1 | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom1      |                                                                   |
+-------------+-------------------------------------------------------------------+
| cra2,cdec2, | Column names in **cross** sample table (tab2).                    |
| cred2,cwei2 | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom2      |                                                                   |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _indic-rppiA:

Input Parameters Dictionary (rppiA)
===================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('rppiA')                            |
+-------------+-------------------------------------------------------------------+
| h0          | Hubble constant [km/s/Mpc]. Default=100.                          |
+-------------+-------------------------------------------------------------------+
| omegam      | Omega matter. Default=0.3                                         |
+-------------+-------------------------------------------------------------------+
| omegal      | Omega lambda. Default=0.7                                         |
+-------------+-------------------------------------------------------------------+
| autogrid    | If ``autogrid=True`` choose the optimum nr. of cells              |            
|             | ``(mxh1,mxh2,mxh3)`` of the skip table (SK). Default=True         |
+-------------+-------------------------------------------------------------------+
| dens        | Custom nr. of particles per SK cell used when ``autogrid=True``.  |
|             | No need to specify unless desired. Default=None                   |
+-------------+-------------------------------------------------------------------+
| mxh1        | Nr. of DEC cells of the SK table. Only relevant if                |
|             | ``autogrid=False``. Default=30                                    |
+-------------+-------------------------------------------------------------------+
| mxh2        | Nr. of RA cells of the SK table. Only relevant if                 |
|             | ``autogrid=False``. Default=180                                   |
+-------------+-------------------------------------------------------------------+
| mxh3        | Nr. of DCOM cells of the SK table. Only relevant if               |
|             | ``autogrid=False``. See :ref:`notemxh3` below. Default=40         |
+-------------+-------------------------------------------------------------------+
| doboot      | If ``True``, calculate bootstrap counts and error bars.           |
|             | Default=False                                                     |
+-------------+-------------------------------------------------------------------+
| nbts        | Nr. of bootstrap samples. Only relevant if ``doboot=True``.       |
|             | Default=50                                                        |
+-------------+-------------------------------------------------------------------+
| bseed       | Fixed seed for boostrap RNG. Always set ``bseed>0`` if running    |
|             | in paralell. Default=1245                                         |
+-------------+-------------------------------------------------------------------+
| wfib        | If ``True`` apply SDSS fiber correction for pairs closer than     |
|             | 55". See *wfiber()* Fortran function. Default=False               |
+-------------+-------------------------------------------------------------------+
| nsepp       | Nr. of projected separation bins. Default=22                      |
+-------------+-------------------------------------------------------------------+
| seppmin     | Minimum projected separation to consider [Mpc/h]. Default=0.01    |
+-------------+-------------------------------------------------------------------+
| dsepp       | Size of projected bins (in dex if log bins). Default=0.15         |
+-------------+-------------------------------------------------------------------+
| logsepp     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
+-------------+-------------------------------------------------------------------+
| nsepv       | Nr. of radial separation bins. Default=1                          |
+-------------+-------------------------------------------------------------------+
| dsepv       | Size of radial bins [Mpc/h]. Default=40.                          |
+-------------+-------------------------------------------------------------------+
| calcdist    | If ``False``, take comov. distances from input tables             |
|             | instead of calculating them. Default=True                         |
+-------------+-------------------------------------------------------------------+
| file        | File name of data sample. Only informative. Default=''            |
+-------------+-------------------------------------------------------------------+
| description | Short description of the run. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| cra,cdec,   | Column names in th sample table (tab).                            |
| cred,cwei,  | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom       |                                                                   |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _indic-rppiC:

Input Parameters Dictionary (rppiC)
===================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('rppiC')                            |
+-------------+-------------------------------------------------------------------+
| h0          | Hubble constant [km/s/Mpc]. Default=100.                          |
+-------------+-------------------------------------------------------------------+
| omegam      | Omega matter. Default=0.3                                         |
+-------------+-------------------------------------------------------------------+
| omegal      | Omega lambda. Default=0.7                                         |
+-------------+-------------------------------------------------------------------+
| autogrid    | If ``autogrid=True`` choose the optimum nr. of cells              |            
|             | ``(mxh1,mxh2,mxh3)`` of the skip table (SK). Default=True         |
+-------------+-------------------------------------------------------------------+
| dens        | Custom nr. of particles per SK cell used when ``autogrid=True``.  |
|             | No need to specify unless desired. Default=None                   |
+-------------+-------------------------------------------------------------------+
| mxh1        | Nr. of DEC cells of the SK table. Only relevant if                |
|             | ``autogrid=False``. Default=30                                    |
+-------------+-------------------------------------------------------------------+
| mxh2        | Nr. of RA cells of the SK table. Only relevant if                 |
|             | ``autogrid=False``. Default=180                                   |
+-------------+-------------------------------------------------------------------+
| mxh3        | Nr. of DCOM cells of the SK table. Only relevant if               |
|             | ``autogrid=False``. See :ref:`notemxh3` below. Default=40         |
+-------------+-------------------------------------------------------------------+
| doboot      | If ``True``, calculate bootstrap counts and error bars.           |
|             | Default=False                                                     |
+-------------+-------------------------------------------------------------------+
| nbts        | Nr. of bootstrap samples. Only relevant if ``doboot=True``.       |
|             | Default=50                                                        |
+-------------+-------------------------------------------------------------------+
| bseed       | Fixed seed for boostrap RNG. Always set ``bseed>0`` if running    |
|             | in paralell. Default=1245                                         |
+-------------+-------------------------------------------------------------------+
| wfib        | If ``True`` apply SDSS fiber correction for pairs closer than     |
|             | 55". See *wfiber()* Fortran function. Default=False               |
+-------------+-------------------------------------------------------------------+
| nsepp       | Nr. of projected separation bins. Default=22                      |
+-------------+-------------------------------------------------------------------+
| seppmin     | Minimum projected distance to consider [Mpc/h]. Default=0.01      |
+-------------+-------------------------------------------------------------------+
| dsepp       | Size of projected bins (in dex if log bins). Default=0.15         |
+-------------+-------------------------------------------------------------------+
| logsepp     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
+-------------+-------------------------------------------------------------------+
| nsepv       | Nr. of radial separation bins. Default=1                          |
+-------------+-------------------------------------------------------------------+
| dsepv       | Size of radial bins [Mpc/h]. Default=40.                          |
+-------------+-------------------------------------------------------------------+
| calcdist    | If ``False``, take comov. distances from input tables             |
|             | instead of calculating them. Default=True                         |
+-------------+-------------------------------------------------------------------+
| file        | File name of **data** sample. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| file1       | File name of **random** sample. Only informative. Default=''      |
+-------------+-------------------------------------------------------------------+
| description | Short description of the run. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| cra,cdec,   | Column names in **data** sample table (tab).                      |
| cred,cwei,  | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom       |                                                                   |
+-------------+-------------------------------------------------------------------+
| cra1,cdec1, | Column names in **random** sample (tab1).                         |
| cred1,cwei1 | Default=('ra' , 'dec', 'z', 'wei', 'dcom')                        |
| cdcom1      |                                                                   |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _notemxh3:

Note on mxh3
============
Due to perfomance reasons, the number of cells in the radial (comoving) distance
actually used to build the skip table is always set as ``mxh3=int((dcmax-dcmin)/rvmax)``.
Hence, the parameter ``mxh3`` supplied at input will be ignored unless it is 
smaller than this optimum value.

Note however that ``mxh3`` is only relevant for the performance of the algorithms.
It is **not** related with the number of radial bins ``nsepv`` where we want to
get output counts.


.. |br| raw:: html

   <br />
