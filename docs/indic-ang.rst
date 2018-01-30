=============
Angular Space
=============

.. _indic-acf:

Input Parameters Dictionary (acf)
=================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('acf')                              |
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
| pxorder     | Data ordering method. See :ref:`pxorder` for details and          |
|             | options. Default='natural'                                        |
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
| nsept       | Nr. of angular bins. Default=22                                   |
+-------------+-------------------------------------------------------------------+
| septmin     | Minimum separation to consider [deg]. Default=0.01                |
+-------------+-------------------------------------------------------------------+
| dsept       | Size of angular bins (in dex if log bins). Default=0.15           |
+-------------+-------------------------------------------------------------------+
| logsept     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
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
| cwei        | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| cra1,cdec1, | Column names in **random** sample table (tab1).                   |
| cwei1       | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| custRAbound | Specify custom RA boundaries for samples that cross the RA=0      |
|             | limit. See :ref:`custRAbound`. Default=None                       |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _indic-accf:

Input Parameters Dictionary (accf)
==================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('accf')                             |
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
| pxorder     | Data ordering method. See :ref:`pxorder` for details and          |
|             | options. Default='natural'                                        |
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
| nsept       | Nr. of angular bins. Default=22                                   |
+-------------+-------------------------------------------------------------------+
| septmin     | Minimum separation to consider [deg]. Default=0.01                |
+-------------+-------------------------------------------------------------------+
| dsept       | Size of angular bins (in dex if log bins). Default=0.15           |
+-------------+-------------------------------------------------------------------+
| logsept     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
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
| cwei        | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| cra1,cdec1, | Column names in **random** sample table (tab1).                   |
| cwei1       | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| cra2,cdec2, | Column names in **cross** sample table (tab2).                    |
| cwei2       | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| custRAbound | Specify custom RA boundaries for samples that cross the RA=0      |
|             | limit. See :ref:`custRAbound`. Default=None                       |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _indic-thA:

Input Parameters Dictionary (thA)
=================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('thA')                              |
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
| pxorder     | Data ordering method. See :ref:`pxorder` for details and          |
|             | options. Default='natural'                                        |
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
| nsept       | Nr. of angular bins. Default=22                                   |
+-------------+-------------------------------------------------------------------+
| septmin     | Minimum separation to consider [deg]. Default=0.01                |
+-------------+-------------------------------------------------------------------+
| dsept       | Size of angular bins (in dex if log bins). Default=0.15           |
+-------------+-------------------------------------------------------------------+
| logsept     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
+-------------+-------------------------------------------------------------------+
| file        | File name of **data** sample. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| description | Short description of the run. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| cra,cdec,   | Column names in **data** sample table (tab).                      |
| cred,cwei   | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| custRAbound | Specify custom RA boundaries for samples that cross the RA=0      |
|             | limit. See :ref:`custRAbound`. Default=None                       |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+


.. _indic-thC:

Input Parameters Dictionary (thC)
=================================

Dictionary with attribute-style access, that stores all input parameters for the
code, plus some useful runtime information during output (see :ref:`outaddrt`)

+-------------+-------------------------------------------------------------------+
| Parameter   | Description                                                       |
+=============+===================================================================+
| kind        | Kind of correlation function ('thC')                              |
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
| nsept       | Nr. of angular bins. Default=22                                   |
+-------------+-------------------------------------------------------------------+
| septmin     | Minimum separation to consider [deg]. Default=0.01                |
+-------------+-------------------------------------------------------------------+
| dsept       | Size of angular bins (in dex if log bins). Default=0.15           |
+-------------+-------------------------------------------------------------------+
| logsept     | If ``True`` use log-spaced bins. Otherwise use linear-spaced      |
|             | bins. Default=True                                                |
+-------------+-------------------------------------------------------------------+
| file        | File name of data sample. Only informative. Default=''            |
+-------------+-------------------------------------------------------------------+
| file1       | File name of random sample. Only informative. Default=''          |
+-------------+-------------------------------------------------------------------+
| description | Short description of the run. Only informative. Default=''        |
+-------------+-------------------------------------------------------------------+
| pxorder     | Sorting method for data. See :ref:`pxorder` for options.          |
|             | Default='natural'                                                 |
+-------------+-------------------------------------------------------------------+
| cra,cdec,   | Column names in **data** sample table (tab).                      |
| cwei        | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| cra1,cdec1, | Column names in **random** sample table (tab1).                   |
| cwei1       | Default=('ra', 'dec', 'wei')                                      |
+-------------+-------------------------------------------------------------------+
| custRAbound | Specify custom RA boundaries for samples that cross the RA=0      |
|             | limit. See :ref:`custRAbound`. Default=None                       |
+-------------+-------------------------------------------------------------------+
| outfn       | Base name for all ouput files (e.g. /home/myuser/redagn)          |
+-------------+-------------------------------------------------------------------+

