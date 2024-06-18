"""
============================================================
GUNDAM :  A Toolkit For Fast Two-Point Correlation Functions
============================================================

@author: Emilio Donoso <edonoso@conicet.gov.ar>
"""
#   Define imports   ==========================================================
import os
import sys
import time
from collections import OrderedDict
from copy import deepcopy
from logging import INFO, FileHandler, StreamHandler, getLogger

import gundam.cflibfor as cff
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Column, Table
from munch import Munch

# ==============================================================================

cps = 0  # Current plot style, as defined in plotcf(). Incremented each time
# that function is called.


# =============================================================================
def set_threads(t):
    """
    Set maximum number of threads to use. When ``t=-1`` defaults to the value
    returned by multiprocessing.cpu_count()

    .. rubric:: Parameters

    t :  integer
        Number of threads desired

    .. rubric:: Returns

    t : integer
        Number of threads to use
    """
    from multiprocessing import cpu_count

    maxt = cpu_count()
    if t <= 0:
        t = maxt
    else:
        t = min(t, maxt)
    return t


# =============================================================================
def comdis_worker(z, cosmo):
    """Worker function for comoving distancesm For more info check the documentation
    of astropy.cosmology module
    """
    return cosmo.comoving_distance(z).value


# =============================================================================
def comdis(z, par, nproc):
    """
    Calculate comoving distances in parallel using multiprocessing and astropy
    cosmology routines

    .. rubric:: Parameters

    z :  array
        Array of redshifts

    par : Munch dictionary
        Input parameters. Used to extract cosmology values

    nproc : integer
        Number of processors to use

    .. rubric:: Returns

    res : array
        Array of comoving distances
    """
    from multiprocessing import Pool

    from astropy.cosmology import LambdaCDM

    # Create cosmology object
    cosmo = LambdaCDM(H0=par.h0, Om0=par.omegam, Ode0=par.omegal)

    pool = Pool(processes=nproc)  # create pool
    zz = np.array_split(z, nproc)  # split input array

    # Apply comdis_worker() to each chunk of input
    res = [pool.apply_async(comdis_worker, [chk, cosmo]) for chk in zz]
    pool.close()
    pool.join()
    res = [chk.get() for chk in res]
    res = np.concatenate(res)
    return res


# =============================================================================
def cfindex(path="./"):
    """
    List file_name + description of all **count** (.cnt) files in a given path.
    Useful to quickly check a dir with dozens of correlation function runs
    and therefore hundreds of files.

    .. rubric:: Parameters

    path : string
        Path to list descriptions of .cnt files inside
    """
    # Get files
    cnt_files = [f for f in os.listdir(path) if f.endswith(".cnt")]
    cnt_files.sort()
    # Get descriptions
    descr = [readcounts(path + f, silent=True).par.description for f in cnt_files]
    # Print left aligned table
    t = Table([cnt_files, descr], names=["File", "Description"])
    t.pprint(align="<", show_dtype=False, max_width=400)


# =============================================================================
def qprint(self):
    """
    Prints a quick, nicely formatted version of Munch dictionaries, such as the
    **counts** ouput dictionary or the **par** input dictionary used by Gundam
    routines. Very useful to quickly explore what is stored inside.

    This routine is added dynamically to the Munch class, so it can be accessed
    as ``any_munch_obj.qprint()``

    Note **counts** output dictionaries can also be explored using
    :func:`gundam.cnttable` to display a  more elaborated tabular view of
    (exclusively) the counts calculated.

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun

        # Explore the many arrays and counts of a typical run
        cred = gun.readcounts('redgalaxies.cnt')
        cred.qprint()

        # Check the parameters used to create the run
        cred.par.qprint()
    """

    def make_custom_sort(orders):
        # Allows to sort a dict in a custom order. Useful to print nicely long
        # dicts such as the par objects used by Gundam functions. Extracted from :
        #   http://stackoverflow.com/questions/12031482/custom-sorting-python-dictionary
        orders = [{k: -i for (i, k) in enumerate(reversed(order), 1)} for order in orders]

        def process(stuff):
            if isinstance(stuff, dict):
                l = [(k, process(v)) for (k, v) in list(stuff.items())]
                keys = set(stuff)
                for order in orders:
                    if keys.issuperset(order):
                        return OrderedDict(sorted(l, key=lambda x: order.get(x[0], 0)))
                return OrderedDict(sorted(l))
            if isinstance(stuff, list):
                return [process(x) for x in stuff]
            return stuff

        return process

    def eprint(v, idn):
        # Print ellipsed version of 1D/2D/3D array
        vs = v.squeeze()
        ndim = vs.ndim
        if ndim == 1:
            txt = eprint1d(vs)
        if ndim == 2:
            txt = eprint2d(vs, idn)
        if ndim == 3:
            txt = eprint3d(vs, idn)
        return txt

    def eprint1d(v):
        # Print ellipsed version of 1D ndarray
        ellip = " ... "
        txt = "[ " + str(v[0]) + ellip + str(v[v.shape[0] // 2]) + ellip + str(v[-1]) + " ]"
        return txt

    def eprint2d(m, idn):
        # Print ellipsed version of 2D ndarray
        ellip = " ... "
        if m.shape[1] >= 3:
            txt1 = "[ " + str(m[0, 0]) + ellip + str(m[0, m.shape[1] // 2]) + ellip + str(m[0, -1]) + " ]"
            txt2 = " " * idn + "[" + ellip + ellip + "]"
            txt3 = (
                " " * idn
                + "[ "
                + str(m[-1, 0])
                + ellip
                + str(m[-1, m.shape[1] // 2])
                + ellip
                + str(m[-1, -1])
                + " ]"
            )
        else:
            txt1 = "[ " + str(m[0, :]) + " ]"
            txt2 = " " * idn + "[" + ellip + "]"
            txt3 = " " * idn + "[ " + str(m[-1, :]) + " ]"
        return "\n".join([txt1, txt2, txt3])

    def eprint3d(m, idn):
        # Print ellipsed version of 3D ndarray
        ellip = " ... "
        txt1 = "[[ " + str(m[0, 0, 0]) + ellip + "///" + ellip + str(m[-1, -1, -1]) + " ]]"
        return "\n".join([txt1])

    spc = " "
    keys = self.keys()

    # Find out if we have a counts or a par object to print accordingly
    if "par" in keys:
        obj_counts = True
        obj_pars = False
    else:
        obj_counts = False
        obj_pars = True

    # We have a counts object
    if obj_counts:
        lj = 12  # nr of characters for left justification of keys
        if self.par.kind == "rppiA":
            msg = "(rp,pi) (auto) counts"
        elif self.par.kind == "rppiC":
            msg = "(rp,pi) (cross) counts"
        elif self.par.kind == "pcf":
            msg = "Projected Correlation"
        elif self.par.kind == "pccf":
            msg = "Projected Cross-correlation"
        elif self.par.kind == "sA":
            msg = "3D (auto) counts"
        elif self.par.kind == "sC":
            msg = "3D (cross) counts"
        elif self.par.kind == "rcf":
            msg = "3D Correlation"
        elif self.par.kind == "rccf":
            msg = "3D Cross-correlation"
        elif self.par.kind == "thA":
            msg = "Angular (auto) counts"
        elif self.par.kind == "thC":
            msg = "Angular (cross) counts"
        elif self.par.kind == "acf":
            msg = "Angular Correlation"
        elif self.par.kind == "accf":
            msg = "Angular Cross-correlation"
        headbar = "=================  " + msg + "  ================="
        print(headbar)
        if self.par.description:
            print("Description".ljust(lj) + "::", self.par.description)

        k = "npt"
        if k in keys:
            print(k.ljust(lj) + "::", self.npt)
        k = "npt1"
        if k in keys:
            print(k.ljust(lj) + "::", self.npt1)
        k = "npt2"
        if k in keys:
            print(k.ljust(lj) + "::", self.npt2)

        allkeys = [
            "rpl",
            "rpm",
            "rpr",
            "sl",
            "sm",
            "sr",
            "thl",
            "thm",
            "thr",
            "wrp",
            "wrperr",
            "xis",
            "xiserr",
            "wth",
            "wtherr",
            "dd",
            "rr",
            "dr",
            "cd",
            "cr",
            "bdd",
            "bcd",
            "rppi",
            "intpi",
            "brppi",
            "intpib",
        ]
        for k in allkeys:
            if k in keys:
                print(k.ljust(lj) + "::", eprint(self[k], lj + 3))

        k = "log"
        if k in keys:
            print(k.ljust(lj) + "::", '"', self.log[0 : self.log.find("\n")], '... "')
        k = "logfortran"
        if k in keys:
            print(k.ljust(lj) + "::", '"', self.logfortran[0 : self.logfortran.find("\n")], '... "')
        k = "par"
        if k in keys:
            print(k.ljust(lj) + ":: {")
            txt = self.par.toJSON(indent=lj + 3, sort_keys=True)
            leng = txt.find("\n", 100)  # just print first few keys
            print(txt[2:leng], "\n" + spc * (lj + 3) + "...", "\n" + spc * (lj + 3) + "}")

        # Print remaining keys, if any
        sprintedkeys = ["npt", "npt1", "npt2", "log", "logfortran", "par"]
        remkeys = set(keys) - set(allkeys) - set(sprintedkeys)
        for k in remkeys:
            v = self[k]
            if type(v) is str:
                v = '"' + v + '"'
            print(k.ljust(lj) + ":: " + str(v))

    # We have a par object
    if obj_pars:
        lj = 15  # nr of characters for left justification of keys
        # This is the order chosen for printing so many parameters
        allkeys = [
            "description",
            "kind",
            "estimator",
            "file",
            "file1",
            "file2",
            "outfn",
            "h0",
            "omegam",
            "omegal",
            "calcdist",
            "autogrid",
            "dens",
            "mxh1",
            "mxh2",
            "mxh3",
            "nsepp",
            "dsepp",
            "seppmin",
            "logsepp",
            "nsepv",
            "dsepv",
            "nseps",
            "dseps",
            "sepsmin",
            "logseps",
            "nsept",
            "dsept",
            "septmin",
            "logsept",
            "doboot",
            "nbts",
            "bseed",
            "wfib",
            "cra",
            "cdec",
            "cred",
            "cwei",
            "cdcom",
            "cra1",
            "cdec1",
            "cred1",
            "cwei1",
            "cdcom1",
            "cra2",
            "cdec2",
            "cred2",
            "cwei2",
            "cdcom2",
            "custRAbound",
        ]

        # Print keys above if present
        for k in allkeys:
            if k in keys:
                v = self[k]
                if type(v) is str:
                    v = '"' + v + '"'
                print(k.ljust(lj) + ":: " + str(v))

        # Print remaining keys, if any
        remkeys = set(keys) - set(allkeys)
        for k in remkeys:
            v = self[k]
            if type(v) is str:
                v = '"' + v + '"'
            print(k.ljust(lj) + ":: " + str(v))

        ## Another way of printing in custom order
        # csort = ['description','kind','estimator','file','file1','outfn','h0',
        #         'omegam','omegal','mxh1','mxh2','mxh3','nsepp','dsepp',
        #         'seppmin','logsepp','nsepv','dsepv',
        #         'doboot','nbts','bseed','wfib','cra','cdec','cred','cwei',
        #         'cra1','cdec1','cred1','cwei1']
        #
        # cs    = make_custom_sort([csort])
        # pards = cs(pard)
        #
        # for k,v in zip(pards.keys(),pards.values()):
        #    if type(v) is str : v = '"' + v + '"'
        #    print(k.ljust(15) + ':: ' + str(v))


# =============================================================================
# Add qprint() method to Munch class. There must be better ways to do this
Munch.qprint = qprint
# =============================================================================


# =============================================================================
def allequal(v):
    """
    Fast way to check if all elements of a 1D-array have the same value. Useful
    to detect when all weights are set to 1.0, and hence to call faster
    versions of the counting routines

    .. rubric:: Parameters

    v :  array_like
        Array to be checked

    .. rubric:: Returns

    res : bool
        True if all elements have the same value
    """
    res = True if (len((v - v[0]).nonzero()[0]) == 0) else False
    return res


# =============================================================================
def addpvcols(x, cfo, basecolname, **kwargs):
    """
    Auxiliary function used by :func:`gundam.cnttable` to append columns to a
    table populated with the fields of **counts** ouput dictionary that store
    pair counts, all with nice column names. Works with 1d counts or 2d counts
    arrays (e.g. those from a pcf run when ``nsepv>1``).

    .. rubric:: Parameters

    x : astropy table
        Table to add data
    cfo : Much dictionary
        Counts dictionary with the count arrays, e.g. cfo.dd, cfo.rr, etc.
    basecolname : string
        The name of the field to add, e.g. `dd`, and also the prefix for the
        column name, which if needed will be appended with `_001`, `_002`, etc.
        for each radial bin
    kwargs :
        Any [key]=value pair to pass to the astropy Column constructor. Intended
        to pass a format specification for the column, such as ``format='.4f'``

    .. rubric:: Returns

    None, it modifies the input table ``x``
    """
    nsepv = cfo.par.nsepv
    for i in range(nsepv):
        if nsepv > 1:
            colname = basecolname + "_" + format(i, "03")
        else:
            colname = basecolname
        colvals = cfo[basecolname][:, i]
        x.add_column(Column(name=colname, data=colvals, **kwargs))


# =============================================================================
def cnttable(cfi, fmt1=None, fmt2=None, write=False, browser=True):
    """
    Shows a nicely formatted tabular view of the count data stored in a
    **counts** output dictionary. The table is printed in `stdout` and optionally
    displayed in the default web browser.

    .. rubric:: Parameters

    cfi : string or Munch dictionary
        Filepath for the counts (.cnt) file, or the **counts** dictionary itself
    fmt1 : string
        Ouput format of numeric fields (bins, correlations and errors). Default='.4f'
    fmt2 : string
        Ouput format of numeric fields (counts). Default='.2f'
    write : bool
        If ``write=True``, save table to disk. Filename will be asked for. Default=False
    browser : bool
        If ``browser=True``, display HTML table in browser. Default=True

    .. rubric:: Returns

    tab : astropy table
        Single table with all relevant counts as columns. Use ``print(tab)`` or
        ``tab.show_in_browser()``

    .. rubric:: Examples

    .. code-block:: python

        # show info for a projected correlation from file on disk
        cnttable('/proj/galred.cnt')
        # show info a from variable in the session's memory
        cnttable(galred)
    """

    if type(cfi) == str:
        cf = readcounts(cfi)  # data comes from file
    else:
        cf = cfi  # data comes from Munch object in memory

    kind = cf.par.kind

    x = Table()  # create table

    fo1 = ".4f" if fmt1 == None else fmt1  # format for bins, cfs and errror
    fo2 = ".2f" if fmt2 == None else fmt2  # format for counts

    if kind in ["pcf", "pccf"]:
        x.add_column(Column(name="rpl", data=cf.rpl, format=fo1))
        x.add_column(Column(name="rpm", data=cf.rpm, format=fo1))
        x.add_column(Column(name="rpr", data=cf.rpr, format=fo1))
        x.add_column(Column(name="wrp", data=cf.wrp, format=fo1))
        x.add_column(Column(name="wrperr", data=cf.wrperr, format=fo1))
        if kind == "pcf":
            addpvcols(x, cf, "dd", format=fo2)
            addpvcols(x, cf, "rr", format=fo2)
            if "dr" in cf:
                addpvcols(x, cf, "dr", format=fo2)
        if kind == "pccf":
            addpvcols(x, cf, "cd", format=fo2)
            addpvcols(x, cf, "cr", format=fo2)

    if kind in ["rcf", "rccf"]:
        x.add_column(Column(name="sl", data=cf.sl, format=fo1))
        x.add_column(Column(name="sm", data=cf.sm, format=fo1))
        x.add_column(Column(name="sr", data=cf.sr, format=fo1))
        x.add_column(Column(name="xis", data=cf.xis, format=fo1))
        x.add_column(Column(name="xiserr", data=cf.xiserr, format=fo1))
        if kind == "rcf":
            x.add_column(Column(name="dd", data=cf.dd, format=fo2))
            x.add_column(Column(name="rr", data=cf.rr, format=fo2))
            if "dr" in cf:
                x.add_column(Column(name="dr", data=cf.dr, format=fo2))
        if kind == "rccf":
            x.add_column(Column(name="cd", data=cf.cd, format=fo2))
            x.add_column(Column(name="cr", data=cf.cr, format=fo2))

    if kind in ["acf", "accf"]:
        x.add_column(Column(name="thl", data=cf.thl, format=fo1))
        x.add_column(Column(name="thm", data=cf.thm, format=fo1))
        x.add_column(Column(name="thr", data=cf.thr, format=fo1))
        x.add_column(Column(name="wth", data=cf.wth, format=fo1))
        x.add_column(Column(name="wtherr", data=cf.wtherr, format=fo1))
        if "dr" in cf:
            x.add_column(Column(name="dr", data=cf.dr, format=fo2))
        if kind == "accf":
            x.add_column(Column(name="cd", data=cf.cd, format=fo2))
            x.add_column(Column(name="cr", data=cf.cr, format=fo2))

    if kind in ["thA", "thC"]:
        x.add_column(Column(name="thl", data=cf.thl, format=fo1))
        x.add_column(Column(name="thm", data=cf.thm, format=fo1))
        x.add_column(Column(name="thr", data=cf.thr, format=fo1))
        if kind == "thA":
            x.add_column(Column(name="dd", data=cf.dd, format=fo2))
        if kind == "thC":
            x.add_column(Column(name="dr", data=cf.dr, format=fo2))

    if kind in ["sA", "sC"]:
        x.add_column(Column(name="sl", data=cf.sl, format=fo1))
        x.add_column(Column(name="sm", data=cf.sm, format=fo1))
        x.add_column(Column(name="sr", data=cf.sr, format=fo1))
        if kind == "sA":
            x.add_column(Column(name="dd", data=cf.dd, format=fo2))
        if kind == "sC":
            x.add_column(Column(name="dr", data=cf.dr, format=fo2))

    if kind in ["rppiA", "rppiC"]:
        x.add_column(Column(name="rpl", data=cf.rpl, format=fo1))
        x.add_column(Column(name="rpm", data=cf.rpm, format=fo1))
        x.add_column(Column(name="rpr", data=cf.rpr, format=fo1))
        if kind == "rppiA":
            addpvcols(x, cf, "dd", format=fo2)
        if kind == "rppiC":
            addpvcols(x, cf, "dr", format=fo2)

    if write:
        fndefault = cf.params.outfn + ".table"
        fn = input("Enter file [" + fndefault + "]: ") or fndefault
        x.write(fn, format="ascii.fixed_width", delimiter="")

    if browser:
        x.show_in_browser()  # show html table in the default web browser
    return x


# =============================================================================
def makebins(nsep, sepmin, dsep, logsep):
    """
    Create arrays of bins in which Gundam will count pairs and estimate
    correlation functions, given the number of bins desired, the minimum bin value and the
    chosen bin width.

    Note units are not needed, but should be interpreted according to the input
    parameters

    .. rubric:: Parameters

    nsep : integer
        Number of bins
    sepmin : float
        Minimum bin location
    dsep : float
        Bin width (dex if ``logsep=1``)
    logsep : bool
        If ``True``, do log-space binning. If ``False``, do linear-space binning

    .. rubric:: Returns

    sep : float array
        Bin locations used by Gundam (left-side + extra bin at right limit)
    sepout : float array
        Left, middle and right-side points of each bin. Useful to plot more easily

    .. rubric:: Examples

    .. code-block:: python

        # Create 18 log bins of size=0.2 dex in redshift and projected space
        seps,sepsout = makebins(18,0.01,0.2,1)
        sepp,seppout = makebins(18,0.01,0.2,1)
        # Create 25 bins of size 2 Mpc in radial space, from 0 to 50 Mpc
        sepv = makebins(25,0.,2.,0)[0]
        # Create instead 1 bin of size 50 Mpc, e.g. to work out the pcf integrated from 0 to 50 Mpc
        sepv = makebins(1,0.,50.,0)[0]
    """
    # Return the limits of bins to do correlation function counts given
    # desired nr. of bins, minimum value and bin width
    sep = np.empty(nsep + 1, dtype=np.float64)
    for i in range(0, nsep + 1):
        if logsep:
            sep[i] = sepmin * 10 ** (i * dsep)
        else:
            sep[i] = sepmin + i * dsep
    # For convenience, also return, leftpoint, midpoint and rightpoint of each bin
    sepout = (sep[0:-1], (sep[0:-1] + sep[1:]) / 2, sep[1:])
    return (sep, sepout)


# =============================================================================
def savecounts(cnt, altname=None):
    """
    Save to disk the **counts** output dictionary returned by the main
    counting routines.

    The default file name is ``cnt.par.outfn`` + `.cnt`, which can be overriden
    as ``altname`` + `.cnt`, if supplied

    .. rubric:: Parameters

    cnt : Munch dictionary
        The **counts** object
    altname : string. Default=None
        If supplied, use an alternative file name instead of ``cnt.par.outfn``
    """

    import pickle

    savename = altname if altname is not None else cnt.par.outfn
    with open(savename + ".cnt", "wb") as f:
        pickle.dump(cnt, f, pickle.HIGHEST_PROTOCOL)

    msg = "> COUNTS object saved in    : " + savename + ".cnt"
    try:
        log = getLogger("cf")
        log.info(msg)
    except:
        print(msg)


# =============================================================================
def readcounts(cfile, silent=False):
    """
    Read from disk the **counts** dictionary generated by the main counting
    routines.

    .. rubric:: Parameters

    cfile : string
        Filepath for the counts (.cnt) file
    silent : bool
        If False, print a status message indicating the file was read. Default=False

    .. rubric:: Returns

    counts : Munch dictionary
        The counts object
    """

    import pickle

    with open(cfile, "rb") as f:
        cnt = pickle.load(f)
        if silent is False:
            print("Counts object read from:", cfile)
    return cnt


# =============================================================================
def plotcf(
    x,
    y,
    yerr,
    fac=1.0,
    write=False,
    figfile=None,
    par=None,
    angunit="deg",
    xlabel=None,
    ylabel=None,
    label=None,
    shift=0.0,
    ploterrbar=True,
    fill=False,
    filtneg=False,
    **kwargs,
):
    """
    Plot a correlation function from arrays of `x`, `y` and `yerr` values.
    Both axes are set to log-space and axes labels are selected automatically
    according to the type of correlation (i.e. given by ``par.kind``)

    .. rubric:: Parameters

    x,y,yerr : float arrays
        x, y coordinates and corresponding errors of the correlation function.
        If ``yerr=0`` or  all elements of yerr are <=0 or ``ploterrbar=False``,
        no errorbar is plotted
    fac : float. Default=1.0
        Multiplication factor for ``y`` and ``yerr``
    write : bool. Default=False
        Save the figure to disk (default format is pdf). See :ref:`Notes <notes-plotcf>`
        to save in other graphic formats
    figfile : string. Default=None
        Specify an alternative file name for the figure. If specified, overrides the
        default which is to take it from ``par.outfn``. Do not add extension.
    par : dictionary of type Munch. Default=None
        Used to pass ``outfn`` name to name saved figures when ``write=True``
    angunit : string. Default='deg'
        * 'arcsec' : set ouput axis in arcsec (``x`` values are unscaled)
        * 'arcmin' : set ouput axis in arcmin (``x`` values are scaled as ``x``/60)
        * 'deg' : set ouput axis in degrees (``x`` values are scaled as ``x``/3600)
    xlabel, ylabel : string. Default=None
        X-axis and Y-axis labels. If supplied, they override the default labels
        deduced from ``par.kind``
    label : string. Default=None
        Label for the legend of the curve. If supplied, it will override the
        default label which is the basename of ``par.outfn``. Note you have
        to issue at least a `plt.legend()` to actually display the legend box
    shift : float. Default=0.0
        Fraction of bin size by which ``x`` values are shifted. Useful to slightly
        separate overlapping curves
    ploterrbar : bool. Default=True
        If ``ploterrbar=True``, plot error bars according to ``yerr``
    fill : bool. Default=False
        If ``fill=True``, plot a filled semi-transparent error region instead of
        the usual error bars
    filtneg : bool. Default=False
        If ``filtneg=True``, filter out points where (y-yerr)<0, i.e. those
        with large errors in a log plot
    kwargs : keyword list
        Any extra [key]=value pairs are passed to the underlying
        :func:`matplotlib.pyplot.plot()` routine, except for ``alpha`` which
        is passed to :func:`pyplot.fill_between()`, ``capsize`` which is passed
        to :func:`pyplot.errorbar()`, and ``figformat`` which is passed to
        :func:`pyplot.savefig()`. Use this to customize colors, linestyles, markers, etc.

    .. _notes-plotcf:

    .. rubric:: Notes

    * Sucessive calls cycle between 4 predefined styles (for color, marker,
      linewidth, etc.) that can be overrriden by passing the corresponding
      [key]=value pairs in ``kwargs``
    * The output graphic format can be changed by passing the ``figformat`` key in
      ``kwargs``, e.g. ``figformat='pdf'``. Any format supported by matplotlib
      is valid.

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun
        c1 = gun.readcounts('redgalPCF.cnt')
        c2 = gun.readcounts('redgalRCF.cnt')

        plt.figure(1)   # Plot w(rp)
        gun.plotcf(c1.rpm,c1.wrp,c1.wrperr,par=c1.par)
        plt.figure(2)   # Plot xi(s)
        gun.plotcf(c2.sm,c2.xis,c2.xiserr,par=c2.par,color='yellow')
    """
    global cps
    killsciform = matplotlib.ticker.ScalarFormatter(useOffset=False)

    # In some cases, CF, counts and errors should be multiplied by factor 2
    y = fac * y
    yerr = fac * yerr

    # Filter negative y-values, if requested. Useful to uncluter both ends
    # of a cf plot, where errors are usually large
    if filtneg:
        pos = (y - yerr) > 0.0
        x = x[pos]
        y = y[pos]
        yerr = yerr[pos]

    # Shift x values by fraction of bin, if desired  ---------
    if shift > 0.0:
        bsiz = np.log10(x[1]) - np.log10(x[0])
        dx = bsiz * shift
        x = 10 ** (np.log10(x) + dx)

    # Get graphic format
    figformat = "pdf" if "figformat" not in kwargs else kwargs.pop("figformat")

    # Define styles and choose label -------------------------
    cls = ["red", "blue", "green", "black"] if "color" not in kwargs else [kwargs.pop("color")] * 4
    mks = ["o", "s", "^", "d"] if "marker" not in kwargs else [kwargs.pop("marker")] * len(cls)
    mkss = [4, 4, 4, 4] if "markersize" not in kwargs else [kwargs.pop("markersize")] * len(cls)
    lin = ["-", "-", "-", "-"] if "linestyle" not in kwargs else [kwargs.pop("linestyle")] * len(cls)
    lwd = [2, 2, 2, 2] if "linewidth" not in kwargs else [kwargs.pop("linewidth")] * len(cls)
    # Get alpha for fills
    alpha = 0.2 if "alpha" not in kwargs else kwargs.pop("alpha")
    # Get capsize for errorbars
    capsize = 3 if "capsize" not in kwargs else kwargs.pop("capsize")

    # Chose label for curve
    if par != None:
        lab = os.path.basename(par.outfn)
    else:
        lab = "data"
    if label != None:
        lab = label

    # Choose x-axis and y-axis labels
    if par != None:
        if par.kind in ["pcf", "pccf", "rppiA", "rppiC"]:
            xtit = r"$r_p \ [h^{-1} Mpc]$" if xlabel == None else xlabel
            ytit = r"$w(r_p)$" if ylabel == None else ylabel
        if par.kind in ["rcf", "rccf", "sA", "sC"]:
            xtit = r"$s \ [h^{-1} Mpc]$" if xlabel == None else xlabel
            ytit = r"$\xi(s) [h^{-1} Mpc]$" if ylabel == None else ylabel
        if par.kind in ["acf", "accf", "thA", "thC"]:
            if angunit == "arcsec":
                xtit = r"$\theta \ [\prime\prime]$" if xlabel == None else xlabel
            if angunit == "arcmin":
                xtit = r"$\theta \ [\prime]$" if xlabel == None else xlabel
                x = x / 60.0
            if angunit == "deg":
                xtit = r"$\theta \ [^{\circ}]$" if xlabel == None else xlabel
                x = x / 3600.0
            ytit = r"$w(\theta)$" if ylabel == None else ylabel
        if par.kind in ["rppiA", "rppiC", "thA", "thC", "sA", "sC"]:
            ytit = r"$counts$" if ylabel == None else ylabel
    else:
        xtit = "x" if xlabel == None else xlabel
        ytit = "y" if ylabel == None else ylabel

    # Plot curve  --------------------------------------------
    plt.plot(
        x,
        y,
        color=cls[cps],
        marker=mks[cps],
        markersize=mkss[cps],
        linestyle=lin[cps],
        linewidth=lwd[cps],
        label=lab,
        **kwargs,
    )

    # Plot error bars  ---------------------------------------
    if ((np.shape(yerr) != ()) or ((y > 0.0).any())) and (ploterrbar):
        if not (fill):
            plt.errorbar(
                x,
                y,
                yerr=yerr,
                fmt="none",
                ecolor=cls[cps],
                elinewidth=lwd[cps],
                capsize=capsize,
                capthick=lwd[cps],
                label=None,
            )
        else:
            plt.fill_between(x, y - yerr, y + yerr, where=(y - yerr) > 0.0, alpha=alpha, color=cls[cps])
            # The where>0 condition avoids points going negative in log scale

    # Others  -------------------------------------------
    plt.xlabel(xtit)
    plt.ylabel(ytit)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(killsciform)  # change labels from sci to plain
    ax.yaxis.set_major_formatter(killsciform)

    if write:
        if par != None:
            pat = figfile if figfile is not None else par.outfn
            kind = "" if figfile is not None else par.kind
        else:
            pat = figfile
            kind = ""
        plt.savefig(pat + "." + kind + "." + figformat, format=figformat)

    if cps < len(cls) - 1:
        cps = cps + 1  # increment current plot style
    else:
        cps = 0

    # Update legend if exists. Useful for adding curves to figures made
    # with comparecf()
    if ax.get_legend():
        ax.legend(frameon=False, fontsize="small")


# =============================================================================
def cntplot(cnt, **kwargs):
    """
    Plot a correlation function from a **counts** output dictionary (either read
    from disk or passed directly). Both axes are set to log-space and axes labels
    are selected automatically according to the type of correlation (i.e. given
    by ``par.kind``)

    This is a wrapper for :func:`gundam.plotcf`, so all of its parameters can
    be specified too.

    .. rubric:: Parameters

    cnt : string or Munch dictionary
        Filepath for the counts (.cnt) file, or the **counts** dictionary itself
    kwargs : keyword list
        Any extra [key]=value pairs are passed to the underlying
        :func:`gundam.plotcf` routine

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun

        # Read a pcf run and plot the correlation function
        cnt1 = gun.readcounts('/p01/redgalsP.cnt')
        cntplot(cnt1)
        # Plot the correlation function from a .cnt file
        cntplot('/p01/redgalsA.cnt', label='angcf of redgals', fill=True)
    """

    # If cnt is a filepath, load data from file
    if type(cnt) == str:
        cnt = readcounts(cnt)

    # Find out the kind of correlation
    kind = cnt.par.kind
    # Choose data keys and titles to plot
    if kind in ["acf", "accf"]:
        x = cnt.thm
        y = cnt.wth
        yerr = cnt.wtherr
    elif kind in ["pcf", "pccf"]:
        x = cnt.rpm
        y = cnt.wrp
        yerr = cnt.wrperr
    elif kind in ["rcf", "rccf"]:
        x = cnt.sm
        y = cnt.xis
        yerr = cnt.xiserr
    elif kind == "rppiA":
        x = cnt.rpm
        y = cnt.intpi
        yerr = 0
    elif kind == "thA":
        x = cnt.thm
        y = cnt.dd
        yerr = 0
    elif kind == "sA":
        x = cnt.sm
        y = cnt.dd
        yerr = 0
    elif kind == "rppiC":
        x = cnt.rpm
        y = cnt.intpi
        yerr = 0
    elif kind == "thC":
        x = cnt.thm
        y = cnt.dr
        yerr = 0
    elif kind == "sC":
        x = cnt.sm
        y = cnt.dr
        yerr = 0

    # Do the plot
    plotcf(x, y, yerr, par=cnt.par, **kwargs)


# =============================================================================
def comparecf(
    clist1,
    clist2=None,
    shift=0.0,
    fac=1.0,
    ploterrbar=True,
    fill=False,
    filtneg=False,
    label=None,
    plotratio=False,
    ratioxrange=None,
    color1=None,
    marker1=None,
    markers1=None,
    linestyle1=None,
    linewidth1=None,
    color2=None,
    marker2=None,
    markers2=None,
    linestyle2=None,
    linewidth2=None,
    f=None,
    ax1=None,
    ax2=None,
):
    """
    Plot multiple correlation functions in a single figure for easy comparison.
    Optionally show an additional lower panel displaying the ratio of each
    function respect to a single "control" correlation (many-to-one) or to
    multiple correlations (one-to-one).

    .. rubric:: Parameters

    clist1 :  list of Munch dictionaries / list of strings
        Count dictionaries or filepaths of .cnt files of correlations
    clist2 : list of Munch dictionaries / list of strings. Default=None
        List of control correlations. When ``plotratio=True``, the y-values
        of each correlation curve in ``clist1`` are divided by those in ``clist2``
        (one-to-one) and plotted in a lower panel. If ``clist2`` has a single
        element, the ratios are from all ``clist1`` divided by the single ``clist1``.
        See :ref:`Notes <notes-comparecf>` for more details
    shift : float. Default=0.0
        Fraction of bin size by which ``x`` values are shifted. Useful to slightly
        separate overlapping curves
    fac : float. Default=1.0
        Multiplication factor for ``y`` and ``yerr``
    ploterrbar : bool. Default=True
        If ``ploterrbar=True``, plot error bars according to ``yerr``
    fill : bool. Default=False
        If ``fill=True``, plot a filled semi-transparent error region instead of
        the usual error bars
    filtneg : bool. Default=False
        If ``filtneg=True``, filter out points where (y-yerr)<0, i.e. those
        with large errors in a log plot
    label: list
        Optional list of strings to label each correlation function. If ommited,
        the values are taken from the ``outfn`` key stored in the count objects
    plotratio : bool. Default=False
        If ``plotratio=True``, plot also a lower panel with ratios of
        ``clist1`` respect to ``clist2``
    ratioxrange : list of form [xmin,xmax]
        Only plot the ratio between xmin and xmax
    color1,marker1,markers1,linestyle1,linewidth1 : lists
        List of colors, marker symbols, marker sizes, line styles and line widths
        for curves in ``clist1``
    color2,marker2,markers2,linestyle2,linewidth2 : lists
        List of colors, marker symbols, marker sizes, line styles and line widths
        for control curves in ``clist2``
    f : figure instance
        Handle of existing Figure instance
    ax1,ax2 : axes instances
        Handles of correlation plot and ratio plot axes

    .. _notes-comparecf:

    .. rubric:: Notes

    The correlation curves in ``clist2`` are **not** plotted in the correlation
    function panel while curves present in **both clists** are **not shown**
    in the ratio panel (i.e. to avoid ratios of curves respect to themselves).

    .. rubric:: Returns

    (f,ax1,ax2) or (f,ax1) : tuple of handles
        Handles of figure, correlation axis (ax1) and ratio axis (ax2), if present

    .. rubric:: Examples

    .. code-block:: python

        # Compare two w(rp) correlations
        comparecf(['galred.cnt', 'galblue.cnt'])
        # Compare one acf on disk with another passed as a counts ouput dictionary
        comparecf(['/proj/galred.cnt', qso])
        # Compare multiple samples and plot sample/control ratios
        f,ax1,ax2 = comparecf(['galred.cnt', 'galblue.cnt', 'galgreen.cnt'], clist2=['allgals.cnt'], fill=True, plotratio=True)
        # Add another curve to previous plot
        comparecf(['qso.cnt'], clist2=['control_qso.cnt'], color2=['k'], f=f, ax1=ax1, ax2=ax2, plotratio=True)
    """
    from matplotlib.pyplot import legend

    n1 = len(clist1)
    n2 = len(clist2) if clist2 else 0

    # Define styles for cfs and ratios  -----------------
    cls1 = color1 if color1 else ["red", "blue", "green", "black"]
    mks1 = marker1 if marker1 else ["o", "s", "^", "d"]
    mkss1 = markers1 if markers1 else [4, 4, 4, 4]
    lin1 = linestyle1 if linestyle1 else ["-", "-", "-", "-"]
    lwd1 = linewidth1 if linewidth1 else [2, 2, 2, 2]

    cls2 = color2 if color2 else cls1
    mks2 = marker2 if marker2 else mks1
    mkss2 = markers2 if markers2 else mkss1
    lin2 = linestyle2 if linestyle2 else lin1
    lwd2 = linewidth2 if linewidth2 else lwd1

    # Set up 1 or 2 axis  ----------------------------
    if plotratio:
        if ax1 is None:
            f, (ax1, ax2) = plt.subplots(
                2, 1, sharex="col", gridspec_kw={"height_ratios": [1.5, 1], "hspace": 0.05}
            )
        kk = f.sca(ax1)  # set current axis to 1st axis

    # Find out the kind of correlation ---------------
    kind = readcounts(clist1[0]).par.kind if type(clist1[0]) == str else clist1[0].par.kind

    for i in range(n1):
        # Data comes from file or Munch object in variable
        data = readcounts(clist1[i]) if type(clist1[i]) == str else clist1[i]
        if kind == "pcf":
            x = data.rpm
            y = data.wrp
            yerr = data.wrperr
            bsiz = data.par.dsepp
        if kind == "xis":
            x = data.mids
            y = data.xis
            yerr = data.xiserr
            bsiz = data.par.dseps
        if kind == "wth":
            x = data.midth
            y = data.wth
            yerr = data.wtherr
            bsiz = data.par.dsep
        if shift > 0:
            dx = bsiz * shift * i
            x = 10 ** (np.log10(x) + dx)
        if label is not None:
            lab = label[i]
        else:
            lab = label
        pars = data.par

        # Plot the ith correlation function  ----------------- kind=kind
        plotcf(
            x,
            y,
            yerr,
            fac=fac,
            write=False,
            par=pars,
            ploterrbar=ploterrbar,
            fill=fill,
            label=lab,
            filtneg=filtneg,
            color=cls1[i],
            marker=mks1[i],
            markersize=mkss1[i],
            linestyle=lin1[i],
            linewidth=lwd1[i],
        )
        legend(frameon=False, fontsize="small")

        # Plot the ratio of corr.func. Single control case  -----
        if plotratio and (n2 == 1):
            pos2in1 = [j for j, k in enumerate(clist1) if k == clist2[0]]
            if pos2in1 != [i]:
                kk = f.sca(ax2)  # set current axis to 2nd axis
                # Data comes from file or Munch object in variable
                data2 = readcounts(clist2[0]) if type(clist2[0]) == str else clist2[0]
                if kind == "pcf":
                    yc = data2.wrp
                if kind == "xis":
                    yc = data2.xis
                if kind == "wth":
                    yc = data2.wth
                # Do ratio --------------------------------------
                yy = y / yc
                # Limit ratio to chosen x-axis range  -----------
                if ratioxrange:
                    good = (x > ratioxrange[0]) & (x < ratioxrange[1])
                    x = x[good]
                    yy = yy[good]
                # Plot ratio  -----------------------------------
                plt.plot(
                    x,
                    yy,
                    color=cls2[i],
                    marker=mks2[i],
                    markersize=mkss2[i],
                    linestyle=lin2[i],
                    linewidth=lwd2[i],
                )
                plt.hlines(
                    1.0, plt.gca().get_xlim()[0], plt.gca().get_xlim()[1], linestyles="dashed", linewidth=1.0
                )
                plt.xlabel(ax1.get_xlabel())  # get x-axis label as defined by plotcf()

                kk = f.sca(ax1)  # go back to 1st axis

        # Plot the ratio of corr.func. Multiple control case  -----
        if plotratio and (n2 > 1):
            raise NameError("Multiple control case not yet simplemented")

    if plotratio:
        ax1.set_xlabel("")  # remove x-label from 1st axis
    ret = (f, ax1, ax2) if plotratio else (plt.gcf(), plt.gca())
    return ret


# =============================================================================
def fitpowerlaw(x, y, yerr, iguess=[1.0, -1.0], fitrange=None, plot=False, markfitrange=False, **kwargs):
    r"""
    Fit a power-law of the form :math:`ax^{\gamma}` to a correlation function over
    a given x-coordinate range. Optionally plot the fitted curve

    .. rubric:: Parameters

    x : float array
        x-coordinates such as cnt.rpm, cnt.thm, etc.
    y : float array
        x-coordinates such as cnt.wrp, cnt.wth, etc.
    yerr : float array
        Errors in y-coordinates such as cnt.wrperr, cnt.wtherr, etc.
    iguess : list of floats. Default=[1., -1.]
        Initial guesses for :math:`a` and :math:`\gamma`
    fitrange : float array of form [xmin,xmax]. Default=None
        Restrict fit to points inside the given interval
    plot : bool. Default=False
        Plot the fitted power-law curve
    markfitrange : bool. Default=False
        Overlay marks for points actually used in the fitting
    kwargs : keyword list
        Any extra [key]=value pairs are passed to :func:`matplolib.pyplot.plot()`
        Use this to customize colors, linestyles, markers, etc.

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun

        c1 = gun.readcounts('galaxies.cnt')
        cntplot(c1)
        gun.fitpowerlaw(c1.rpm, c1.wrp, c1.wrperr, plot=True)
    """
    import scipy.optimize

    # Filter out points outside fitrange
    if fitrange:
        idx = np.where((x >= fitrange[0]) & (x <= fitrange[1]))
        x = x[idx]
        y = y[idx]
        yerr = yerr[idx]

    # Take logs
    xx = np.log10(x)
    yy = np.log10(y)
    yyerr = yerr / y

    # Filter out points with nan/infs/etc.
    idx = np.isfinite(yy) & np.isfinite(yyerr)
    xx = xx[idx]
    yy = yy[idx]
    yyerr = yyerr[idx]

    # Define the (line) fitting function
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    # Carry out fitting
    # pinit = [1.0, -1.0]   # initial guess
    out = scipy.optimize.leastsq(errfunc, iguess, args=(xx, yy, yyerr), full_output=1)
    pfinal = out[0]  # final parameters
    covar = out[1]  # covariance matrix

    # Extract fit parameters and errors
    gamma = pfinal[1]
    a = 10.0 ** pfinal[0]
    gammaerr = np.sqrt(covar[0][0])
    aerr = np.sqrt(covar[1][1]) * a

    print("Fitting Power Law --> y = a*x^gamma")
    print("===================================")
    print("  gamma =", gamma, "\u00B1", gammaerr)
    print("  a     =", a, "\u00B1", aerr)
    print("===================================")
    if plot:
        cls = "black" if "color" not in kwargs else kwargs.pop("color")
        lin = "--" if "linestyle" not in kwargs else kwargs.pop("linestyle")
        lwd = 1 if "linewidth" not in kwargs else [kwargs.pop("linewidth")]
        xdom = np.linspace(plt.xlim()[0], plt.xlim()[1], 10)
        ydom = a * xdom**gamma
        plt.plot(xdom, ydom, color=cls, linestyle=lin, linewidth=lwd, zorder=1, **kwargs)
        if markfitrange:
            plt.scatter(10**xx, 10**yy, s=80, facecolors="none", edgecolors="blue")

    return (gamma, gammaerr, a, aerr)


# =============================================================================
def cntplot2D(
    cnt, estimator=None, slevel=5, write=False, figfile=None, xlabel=None, ylabel=None, cmap="jet", **kwargs
):
    r"""
    Plot the 2D correlation function in the projected-radial space
    (:math:`r_p` vs :math:`\pi` space) with optional gaussian smoothing and
    contour levels

    .. rubric:: Parameters

    cnt : string or Munch dictionary
        Filepath for the counts (.cnt) file or the **counts** output dictionary
    estimator : string. Default=None
        Estimator for the correlation function. Any of ('NAT','LS','HAM','DP').
        If ``estimator=None``, then it is taken from ``cnt.par.estimator``
    slevel : float. Default=5
        Smoothing level (namely the size of the Gaussian smothing kernel)
    write : bool. Default=False
        Save the figure to disk (default format is pdf). See :ref:`Notes <notes-cntplot2D>`
        to save in other graphic formats
    figfile : string. Default=None
        Specify an alternative file name for the figure. If ``None``, then choose
        ``cnt.par.outfn`` as default. Do not add extension.
    xlabel, ylabel : string. Default=None
        X-axis and Y-axis labels. If supplied, they override the default labels
        (:math:`r_p \ [h^{-1} Mpc]` and :math:`\pi \ [h^{-1} Mpc]`)
    cmap : string. Default='jet'
        Colormap for the plot
    kwargs : keyword list
        Any extra [key]=value pairs are passed to :func:`matplolib.pyplot.pcolor()`
        Use this to customize shading, edges, alpha, etc.

    .. _notes-cntplot2D:

    .. rubric:: Notes

    * The graphic format can be changed by passing the ``figformat`` key in
      ``kwargs``, e.g. ``figformat='pdf'``. Any format supported by matplotlib
      is valid.


    .. rubric:: Examples

    .. code-block:: python

        # Check some nice Fingers of God and the Kaiser squashing
        cntplot2D('lum_red_gals.cnt', cmap='viridis')
    """

    def gaussKern(size):
        """
        Calculate a normalised Gaussian kernel to apply as a smoothing function

        Parameters
        size (int) : kernel size (how many points will be used in the smoothing operation)

        Returns
        g (array(size,size)) : normalised 2D kernel array for use in convolutions
        """
        size = int(size)
        x, y = np.mgrid[-size : size + 1, -size : size + 1]
        g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
        return g / g.sum()

    def smooth(im, n=15):
        """
        Smooth a 2D array im by convolving with a Gaussian kernel of size n

        Parameters
        im, n : 2D array, kernel size

        Returns
        improc : smoothed array (same dimensions as the input array)
        """
        from scipy import signal

        g = gaussKern(n)
        improc = signal.convolve2d(im, g, mode="same", boundary="symm")
        return improc

    # If cnt is a filepath, load data from file
    if type(cnt) == str:
        cnt = readcounts(cnt)

    # Get graphic format
    figformat = "pdf" if "figformat" not in kwargs else kwargs.pop("figformat")

    # Build 2D matrix of correlation function
    est = estimator if estimator is not None else cnt.par.estimator
    nd, nr = cnt.npt * 1.0, cnt.npt1 * 1.0  # for easy typing
    if est == "NAT":
        nf = (nr / nd) * ((nr - 1) / (nd - 1))  # normalizing factor
        xi = nf * (1.0 * cnt.dd / cnt.rr) - 1
    elif est == "LS":
        nf1 = (nr / nd) * ((nr - 1) / (nd - 1))  # normalizing factor 1
        nf2 = (nr - 1) / (2.0 * nd)  # normalizing factor 2
        xi = (nf1 * cnt.dd - nf2 * 2.0 * cnt.dr + cnt.rr) / cnt.rr
    elif est == "HAM":
        nf = 4 * (nd / (nd - 1)) * (nr / (nr - 1))  # normalizing factor
        xi = nf * cnt.dd * cnt.rr / cnt.dr**2 - 1
    elif est == "DP":
        nf = (2.0 * nr) / (nd - 1)  # normalizing factor
        xi = nf * cnt.dd / cnt.dr - 1

    # Take log, replicate over four quadrants and smooth
    logxi = np.log10(xi)
    qTR = logxi.copy()
    qTR[~np.isfinite(qTR)] = 0.0
    qTR = qTR.T  # top right quadrant
    qTL = np.fliplr(qTR)  # top left quadrant
    qBL = np.flipud(qTL)  # bottom left quadrant
    qBR = np.fliplr(qBL)  # bottom right quadrant
    qT = np.hstack((qTL, qTR))  # top half
    qB = np.hstack((qBL, qBR))  # bottom half
    qq = np.vstack((qB, qT))  # full array
    qqs = smooth(qq, n=slevel)  # smoothed full array

    # Extend bins
    sepp = makebins(cnt.par.nsepp, cnt.par.seppmin, cnt.par.dsepp, cnt.par.logsepp)[0]
    sepv = makebins(cnt.par.nsepv, 0.0, cnt.par.dsepv, False)[0]
    exsepv = np.concatenate([-1 * sepv[cnt.par.nsepv : 0 : -1], sepv])
    exsepp = np.concatenate([-1 * sepp[cnt.par.nsepp : 0 : -1], sepp])
    # Get bin coordinates
    x, y = np.meshgrid(exsepp, exsepv)
    # Plot array
    plt.pcolor(x, y, qqs, cmap=cmap, **kwargs)
    # Plot contours
    lev = np.linspace(np.amin(qqs), np.amax(qqs), 15)
    plt.contour(x[0:-1, 0:-1], y[0:-1, 0:-1], qqs, levels=lev, colors="k", linestyles="solid", linewidths=1)
    # Plot titles
    xtit = r"$r_p \ [h^{-1} Mpc]$" if xlabel is None else xlabel
    ytit = r"$\pi \ [h^{-1} Mpc]$" if ylabel is None else ylabel
    plt.xlabel(xtit)
    plt.ylabel(ytit)

    if write:
        pat = figfile if figfile is not None else cnt.par.outfn
        plt.savefig(pat + ".2DCF." + figformat, format=figformat)


# =============================================================================
def writeasc_cf(lb, mb, rb, f, ferr, par, fmt="%17.5f", altname=None):
    """
    Write an ASCII file for w(rp) / xi(s) / w(th) ouput counts produced by the code

    .. rubric:: Parameters

    lb,mb,rb : float arrays
        Left, mid and rigth-side of bins
    f, ferr : float arrays
        Correlation function and its error
    par : Munch dictionary
        Used to pass ``par.outfn`` to name the output file
    fmt : string. Default='%17.5f'
        Numeric formatting string
    altname : string. Default=None
        If supplied, use an alternative file name instead of ``par.outfn``

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun
        c1 = gun.readcounts('redgals.cnt')

        writeasc_cf(c1.rpl, c1.rpm, c1.rpr, c1.wrp, c1.wrperr, c1.par)
    """
    savename = altname if altname is not None else par.outfn

    log = getLogger("cf")
    if par.kind in ["pcf", "pccf"]:
        extension = ".wrp"
        kopf = "         lb              mb              rb            wrp               errwrp"
        savename = savename + extension if altname is None else altname
        msg = "> w(rp) PCF saved in        : " + savename
    elif par.kind in ["rcf", "rccf"]:
        extension = ".xis"
        kopf = "         lb              mb              rb            xis               xiserr"
        savename = savename + extension if altname is None else altname
        msg = "> xi(s) RSCF saved in       : " + savename
    elif par.kind in ["acf", "accf"]:
        extension = ".wth"
        kopf = "         lb              mb              rb            wth               errwth"
        savename = savename + extension if altname is None else altname
        msg = "> w(th) ACF saved in        : " + savename
    m = np.array([lb, mb, rb, f, ferr])
    m = m.transpose()
    np.savetxt(savename, m, fmt=fmt, header=kopf)
    log.info(msg)


# =============================================================================
def writeasc_rppicounts(lb, mb, rb, rppi, par, fmt="%17.5f", cntid=None, altname=None):
    """
    Write an ASCII file for a rp-pi count array produced by the code. This is
    a 2D array of counts in the projected (rp) and radial (pi) directions.

    The columns in the output will be [`lb mb rb tot_counts rppi`] where the
    first 3 are the left, mid and right-side of bins, `tot_counts` are
    the counts integrated for all radial bins, and `rppi` has one column for
    each radial bin

    .. rubric:: Parameters

    lb,mb,rb : float
        Left, mid and rigth-side of bins
    rppi : float array
        2-dimensional count array. Usually this is one of the fields cnt.dd,
        cnt.rr, etc. of a projected correlation run
    par : Munch dictionary
        Used to pass various data, including ``par.outfn`` to name the output file
    fmt : string. Default='%17.5f'
        Numeric formatting string
    cntid : string. Default=None
        ID string for column headers. Usually can be 'dd', 'rr', 'dr', etc.
        Also appended as the extension of the ouput file (when ``altname=None``)
    altname : string. Default=None
        If supplied, use an alternative file name instead of ``par.outfn`` + `.cntid`

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun
        c1 = gun.readcounts('redgals.cnt')

        c1.par.outfn
        '/home/myuser/sdss/redgals'

        # Write the DD counts in rp-pi dimensions
        gun.writeasc_rppicounts(c1.rpl, c1.rpm, c1.rpr, c1.dd, c1.par, cntid='dd')

        # Inspect the output file
        with open('redgals.dd', 'r') as f:
            print(f.read(), end="")

        #  lb        mb        rb        dd         dd_001    dd_002    dd_003    ...
        #  0.10000   0.12417   0.14835   11509.00   2082.00   1500.00   1168.00   ...
        #  0.14835   0.18421   0.22007   20273.00   3122.00   2378.00   1899.00   ...
        #  0.22007   0.27327   0.32647   36169.00   4940.00   3845.00   3283.00   ...
        #  0.32647   0.40539   0.48431   64866.00   8453.00   6302.00   5236.00   ...
        #  ...
    """
    log = getLogger("cf")
    nsepv, nsepp = par.nsepv, par.nsepp
    filename = altname if altname is not None else par.outfn + "." + cntid
    kopf = "         lb              mb               rb            " + cntid
    if nsepv > 1:
        for i in range(1, nsepv + 1):
            kopf = kopf + "          " + cntid + "_" + format(i, "03")
        totc = np.sum(rppi, axis=1)
        totc.shape = (nsepp, 1)

    m = np.transpose(np.array([lb, mb, rb]))
    if nsepv > 1:
        m = np.hstack((m, totc, rppi))
    else:
        m = np.hstack((m, rppi))

    np.savetxt(filename, m, fmt=fmt, header=kopf)
    msg = "> ASCII counts saved in     : " + filename
    log.info(msg)


# =============================================================================
def writeasc_counts(lb, mb, rb, c, par, fmt="%17.5f", cntid=None, altname=None):
    """
    Write an ASCII file for a (1-dimensional) counts array produced by the code.

    The columns in the output will be [`lb mb rb c`] where the
    first 3 are the left, mid and right-side of bins and `c` are the counts

    .. rubric:: Parameters

    lb,mb,rb : float
        Left, mid and rigth-side of bins
    c : float array
        Counts array. Usually this is one of the fields cnt.dd, cnt.dr, cnt.rr,
        etc. of a correlation run
    par : Munch dictionary
        Used to pass ``par.outfn`` to name the output file
    fmt : string. Default='%17.5f'
        Numeric formatting string
    cntid : string. Default=None
        ID string for column header. Usually can be 'dd', 'rr', 'dr', etc.
        Also appended as the extension of the ouput file (when ``altname=None``)
    altname : string. Default=None
        If supplied, use an alternative file name instead of ``par.outfn`` + `.cntid`

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun
        c1 = gun.readcounts('bluegals.cnt')

        # Write the DD counts in angular dimensions. Use an alternative file name
        gun.writeasc_counts(c1.thl, c1.thm, c1.thr, c1.dd, c1.par, cntid='dd', altname='akounts')

        # Inspect the output file
        with open('akounts', 'r') as f:
            print(f.read(), end="")

        # lb        mb        rb           dd
        # 0.01000   0.01206   0.01413     3178.00
        # 0.01413   0.01704   0.01995     6198.00
        # 0.01995   0.02407   0.02818    12765.00
        # 0.02818   0.03400   0.03981    24888.00
        # 0.03981   0.04802   0.05623    49863.00
        # 0.05623   0.06783   0.07943    98883.00
        ...
    """
    log = getLogger("cf")
    filename = altname if altname is not None else par.outfn + "." + cntid
    kopf = "         lb               mb               rb             " + cntid

    m = np.array([lb, mb, rb, c])
    m = m.transpose()
    np.savetxt(filename, m, fmt=fmt, header=kopf)
    msg = "> " + cntid + " counts saved in        : " + filename
    log.info(msg)


# =============================================================================
def savepars(par, altname=None):
    """
    Save the parameters dictionary **par**, such as the one
    generated by :func:`gundam.packpars`, in a JSON file. By default it is
    named as ``par.outfn`` + `.par`

    .. rubric:: Parameters

    par : Munch dictionary
        Input parameters dictionary for Gundam routines
    altname : string. Default=None
        If supplied, use an alternative file name instead of ``par.outfn`` + `.par`

    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun

        # Get default values for an angular CF run and save to disk
        par = gun.packpars(kind='acf', outfn='/proj/acfrun01')
        gun.savepars(par)
    """
    pj = par.toJSON(indent=0)
    filename = altname if altname is not None else par.outfn + ".par"

    with open(filename, "w") as f:
        f.write(pj)
    msg = "> PARAMS saved in           : " + filename
    log = getLogger("cf")
    log.info(msg)


# =============================================================================
def readpars(filename):
    """
    Load from a JSON (.par) file the input parameters dictionary used by many
    Gundam routines.

    .. rubric:: Parameters

    filename : string
        Filepath of .par file
    """
    import munch

    with open(filename, "r") as f:
        s = f.read()
    a = munch.json.loads(s)  # convert string to dict
    b = munch.munchify(a)  # convert dict to bunch
    print("Params read from :", filename)
    return b


# =============================================================================
def tpcf(npt, nrpt, dd, bdd, rr, dr, estimator):
    """
    Return the (auto)correlation function for a given estimator and arrays of
    data and random counts

    If DR counts are not needed (e.g. the 'NAT' estimator), just set ``dr=0``

    If boostrap errors are not needed or available, just set ``bdd`` to a
    zero-valued array with null 2nd dimension, e.g. ``bdd=np.zeros([len(dd),0])``


    .. rubric :: Parameters

    npt,nrpt : integer
        Number of data and random particles
    dd : float array
        DD counts
    bdd : float array
        Bootstrap DD counts
    rr : float array
        RR counts in projected and radial separations
    dr : float array
        DR counts in projected and radial separations
    estimator : string
        Statistical estimator of the correlation function. Default=`NAT`

            * 'NAT' : Natural -> :math:`DD/RR - 1`
            * 'HAM' : Hamilton -> :math:`DD*RR/DR^{2} - 1`
            * 'LS'  : Landy-Szalay -> :math:`(DD - 2DR + RR) / RR`
            * 'DP'  : Davis-Peebles -> :math:`DD/DR - 1`

    .. rubric :: Returns

    xi : float array
        Correlation function
    xierr : float array
        Boostrap error estimate. Set to zero if ``bdd`` is nulled as explained
        above


    .. rubric :: Notes

    See `this paper <http://arxiv.org/pdf/1211.6211v2.pdf>`_ for a nice review on
    estimators and their normalization factors. Here, the normalization factors
    are derived to : (1) keep estimator formulae clean, (2) avoid having
    operations such as (npt*(npt-1)) * dd, where counts are multiplied/divided
    by very big numbers when npt is large.


    .. rubric :: Examples

    .. code-block:: python

        # Calculate the angular CF using the Landy-Szalay estimator
        acf, acferr = gun.tpcf(npt,nrpt,dd,bdd,rr,dr,estimator='LS')
    """
    fn = sys._getframe().f_code.co_name  # get function self name
    msg = "Computing 3d/angular correlation function.....[" + fn + "()]"
    log = getLogger("cf")
    log.info(msg)

    npt, nrpt = float(npt), float(nrpt)
    nseps, nbts = bdd.shape

    # Initialize ouput arrays  -------------------------------------
    xi = np.zeros(nseps)  # the correlation function (initizalize to -1 ?)
    xierr = np.zeros(nseps)  # the error
    bxi = np.zeros(nbts)

    # Loop over bins  ----------------------------------------------
    if estimator == "HAM":
        nf = 4.0 * (npt / (npt - 1)) * (nrpt / (nrpt - 1))  # normalizing factor
        for i in range(nseps):
            if dr[i] > 0.0:
                xi[i] = nf * dd[i] * rr[i] / dr[i] ** 2 - 1.0
                for j in range(nbts):
                    bxi[j] = nf * bdd[i, j] * rr[i] / dr[i] ** 2 - 1.0
                xierr[i] = np.std(bxi) if nbts > 0 else 0.0

    if estimator == "NAT":
        nf = (nrpt / npt) * ((nrpt - 1.0) / (npt - 1.0))  # normalizing factor
        for i in range(nseps):
            if rr[i] > 0.0:
                xi[i] = nf * dd[i] / rr[i] - 1.0
                for j in range(nbts):
                    bxi[j] = nf * bdd[i, j] / rr[i] - 1.0
                xierr[i] = np.std(bxi) if nbts > 0 else 0.0

    if estimator == "LS":
        nf1 = (nrpt / npt) * ((nrpt - 1) / (npt - 1))  # normalizing factor 1
        nf2 = (nrpt - 1.0) / (2.0 * npt)  # normalizing factor 2
        for i in range(nseps):
            if rr[i] > 0.0:
                xi[i] = (nf1 * dd[i] - nf2 * 2.0 * dr[i] + rr[i]) / rr[i]
                for j in range(nbts):
                    bxi[j] = (nf1 * bdd[i, j] - nf2 * 2.0 * dr[i] + rr[i]) / rr[i]
                xierr[i] = np.std(bxi) if nbts > 0 else 0.0

    if estimator == "DP":
        nf = (2.0 * nrpt) / (npt - 1)  # normalizing factor
        for i in range(nseps):
            if dr[i] > 0.0:
                xi[i] = nf * dd[i] / dr[i] - 1.0
                for j in range(nbts):
                    bxi[j] = nf * bdd[i, j] / dr[i] - 1.0
                xierr[i] = np.std(bxi) if nbts > 0 else 0.0

    return (xi, xierr)


# =============================================================================
def tpccf(npt, nrpt, cd, bcd, cr, estimator):
    """
    Return the (cross)correlation function for a given estimator and count
    arrays for data (D), random (R) and cross (C) samples.

    For the moment the only estimator implemented is the Davis-Peebles :
    :math:`\\xi=CD/CR-1`

    If bootstrap errors are not needed or available, just set ``bdd`` to a
    zero-valued array, e.g. ``bdd=np.zeros([len(dd),0])``

    .. rubric:: Parameters

    npt,nrpt : integer
        Number of particles in data (D) and random (R) samples
    cd : float array
        CD counts
    bcd : float array
        Bootstrap CD counts
    cr : float array
        CR counts
    estimator : string
        * 'DP' : Davis-Peebles -> :math:`CD/CR-1`

    .. rubric:: Notes

    C and D are data samples while R is random sample corresponding to D


    .. rubric:: Returns

    fxi : float array
        Cross-correlation function
    fxierr : float array
        Boostrap error estimates


    .. rubric:: Examples

    .. code-block:: python

        import gundam as gun
        c = gun.readcounts('qso_gal.cnt')

        (ccf,ccferr) = tpccf(c.npt, c.nrpt, c.cd, c.bcd, c.cr, estimator='DP')
    """

    fn = sys._getframe().f_code.co_name  # get function self name
    msg = "Computing 3d/angular cross-correlation function....[" + fn + "()]"
    log = getLogger("cf")
    log.info(msg)

    npt, nrpt = float(npt), float(nrpt)
    nseps, nbts = bcd.shape

    # Initialize output arrays  ------------------------------------
    fxi = np.zeros(nseps)
    fxierr = np.zeros(nseps)
    bxi = np.zeros(nbts)

    # Loop over bins  ----------------------------------------------
    if estimator == "DP":
        sf = nrpt / npt  # normalizing factor
        for i in range(nseps):
            if cr[i] > 0.0:
                xi = (cd[i] / cr[i]) * sf - 1.0
                for j in range(nbts):
                    bxi[j] = (bcd[i, j] / cr[i]) * sf - 1.0
                xierr = np.std(bxi) if nbts > 0 else 0.0
            else:
                xi = 0.0  # -1.
                xierr = 0.0  # -1.
            fxi[i] = xi
            fxierr[i] = xierr
    return (fxi, fxierr)


# =============================================================================
def tpccf_wrp(npt, nrpt, cd, bcd, cr, dsepv, estimator):
    """
    Return the projected (cross)correlation function for a given estimator and
    count arrays

    If boostrap errors are not needed or available, just set ``bcd`` to a
    zero-valued array with null 2nd dimension, e.g. ``bcd=np.zeros([len(ddXXX),0])``


    .. rubric :: Parameters

    npt : integer
        Number of data particles
    nrpt : integer
        Number of random particles
    cd : float array
        CD counts in projected and radial separations
    bcd : float array
        Bootstrap CD counts in projected and radial separations
    cr : float array
        CR counts in projected and radial separations
    dsepv : float
        Radial bin size
    estimator : string
        Statistical estimator of the correlation function

            * 'DP'  : Davis-Peebles -> :math:`CD/CR - 1` (C,D data samples, R is random of D)


    .. rubric :: Returns

    wrp : float array
        Projected cross-correlation function
    wrperr : float array
        Boostrap error estimate


    .. rubric :: Examples

    .. code-block:: python

        # remains to do XXXXX
        (wrp,wrperr) = tpccf_wrp(npt,nrpt,cd,bcd,cr,dsepv,estimator='DP')
    """
    fn = sys._getframe().f_code.co_name  # get function self name
    msg = "Computing proj. cross-correlation function.....[" + fn + "()]"
    log = getLogger("cf")
    log.info(msg)

    npt, nrpt = float(npt), float(nrpt)
    nsepp, nsepv = cd.shape
    nbts = bcd.shape[2]

    # Initialize ouput arrays  ------------------------------------
    wrp = np.zeros(nsepp)
    bxi = np.zeros(nbts)
    wrperr = np.zeros(nsepp)

    # Loop over bins  ----------------------------------------------
    if estimator == "DP":
        sf = nrpt / npt  # normalizing factor
        for i in range(nsepp):
            twrp = 0
            tbwrp = 0
            for j in range(nsepv):
                if cr[i, j] > 0.0:
                    xi = (cd[i, j] / cr[i, j]) * sf - 1
                    for k in range(nbts):
                        bxi[k] = (bcd[i, j, k] / cr[i, j]) * sf - 1
                else:
                    xi = 0.0  # -1.
                    bxi = np.zeros(nbts)  # - 1.
                twrp = twrp + xi * dsepv
                tbwrp = tbwrp + bxi * dsepv
            wrp[i] = 2.0 * twrp
            wrperr[i] = 2.0 * np.std(tbwrp) if nbts > 0 else 0.0

    return (wrp, wrperr)


# =============================================================================
def tpcf_wrp(npt, nrpt, dd, bdd, rr, dr, dsepv, estimator):
    """
    Return the projected (auto)correlation function for a given estimator and
    arrays of data and random counts

    If DR counts are not needed (e.g. the 'NAT' estimator), just set ``dr=0``

    If boostrap errors are not needed or available, just set ``bdd`` to a
    zero-valued array with null 2nd dimension, e.g. ``bdd=np.zeros([len(dd),0])``


    .. rubric :: Parameters

    npt,nrpt : integer
        Number of data and random particles
    dd : float array
        DD counts in projected and radial separations
    bdd : float array
        Bootstrap DD counts in projected and radial separations
    rr : float array
        RR counts in projected and radial separations
    dr : float array
        DR counts in projected and radial separations
    dsepv : float
        Bin size in radial direction
    estimator : string
        Statistical estimator of the correlation function

            * 'NAT' : Natural -> :math:`DD/RR - 1`
            * 'HAM' : Hamilton -> :math:`DD*RR/DR^{2} - 1`
            * 'LS'  : Landy-Szalay -> :math:`(DD - 2DR + RR) / RR`
            * 'DP'  : Davis-Peebles -> :math:`DD/DR - 1`

    .. rubric :: Returns

    wrp : float array
        Correlation function
    wrperr : float array
        Boostrap error estimate. Set to zero if ``bdd`` is nulled as explained
        above


    .. rubric :: Notes

    See `this paper <http://arxiv.org/pdf/1211.6211v2.pdf>`_ for a nice review on
    estimators and their normalization factors. Here, the normalization factors
    are derived to : (1) keep estimator formulae clean, (2) avoid having
    operations such as (npt*(npt-1)) * dd, where counts are multiplied/divided
    by very big numbers when npt is large.


    .. rubric:: Examples (xxxx TODO)

    .. code-block:: python

        (wrp,wrperr) = tpcf_wrp(npt,nrpt,ddpv,bddpv,rrpv,drpv,dsepv,estimator='HAM')
    """
    fn = sys._getframe().f_code.co_name  # get function self name
    msg = "Computing proj. correlation function.....[" + fn + "()]"
    log = getLogger("cf")
    log.info(msg)

    npt, nrpt = float(npt), float(nrpt)
    nsepp, nsepv = dd.shape
    nbts = bdd.shape[2]

    # Initialize output arrays  -------------------------------------
    wrp = np.zeros(nsepp)  # the correlation function
    wrperr = np.zeros(nsepp)  # the boostrap error
    bxi = np.zeros(nbts)

    # Loop over bins  ----------------------------------------------
    if estimator == "HAM":
        nf = 4.0 * (npt / (npt - 1.0)) * (nrpt / (nrpt - 1.0))  # normalizing factor
        for i in range(nsepp):
            twrp = 0
            tbwrp = 0
            for j in range(nsepv):
                if dr[i, j] > 0.0:
                    xi = nf * dd[i, j] * rr[i, j] / dr[i, j] ** 2 - 1
                    for k in range(nbts):
                        bxi[k] = nf * bdd[i, j, k] * rr[i, j] / dr[i, j] ** 2 - 1
                else:
                    xi = 0.0  # -1.
                    bxi = np.zeros(nbts)  # - 1.
                twrp = twrp + xi * dsepv
                tbwrp = tbwrp + bxi * dsepv
            wrp[i] = 2.0 * twrp
            wrperr[i] = 2.0 * np.std(tbwrp) if nbts > 0 else 0.0

    if estimator == "NAT":
        nf = (nrpt / npt) * ((nrpt - 1) / (npt - 1))  # normalizing factor
        for i in range(nsepp):
            twrp = 0
            tbwrp = 0
            for j in range(nsepv):
                if rr[i, j] > 0.0:
                    xi = nf * dd[i, j] / rr[i, j] - 1
                    for k in range(nbts):
                        bxi[k] = nf * bdd[i, j, k] / rr[i, j] - 1
                else:
                    xi = 0.0  # -1.
                    bxi = np.zeros(nbts)  # - 1.
                twrp = twrp + xi * dsepv
                tbwrp = tbwrp + bxi * dsepv
            wrp[i] = 2 * twrp
            wrperr[i] = 2 * np.std(tbwrp) if nbts > 0 else 0.0

    if estimator == "LS":
        nf1 = (nrpt / npt) * ((nrpt - 1) / (npt - 1))  # normalizing factor 1
        nf2 = (nrpt - 1) / (2 * npt)  # normalizing factor 2
        for i in range(nsepp):
            twrp = 0
            tbwrp = 0
            for j in range(nsepv):
                if rr[i, j] > 0.0:
                    xi = (nf1 * dd[i, j] - nf2 * 2 * dr[i, j] + rr[i, j]) / rr[i, j]
                    for k in range(nbts):
                        bxi[k] = (nf1 * bdd[i, j, k] - nf2 * 2 * dr[i, j] + rr[i, j]) / rr[i, j]
                else:
                    xi = 0.0  # -1.
                    bxi = np.zeros(nbts)  # - 1.
                twrp = twrp + xi * dsepv
                tbwrp = tbwrp + bxi * dsepv
            wrp[i] = 2.0 * twrp
            wrperr[i] = 2.0 * np.std(tbwrp) if nbts > 0 else 0.0

    if estimator == "DP":
        nf = (2.0 * nrpt) / (npt - 1)  # normalizing factor
        for i in range(nsepp):
            twrp = 0
            tbwrp = 0
            for j in range(nsepv):
                if dr[i, j] > 0.0:
                    xi = nf * dd[i, j] / dr[i, j] - 1
                    for k in range(nbts):
                        bxi[k] = nf * bdd[i, j, k] / dr[i, j] - 1
                else:
                    xi = 0.0  # -1.
                    bxi = np.zeros(nbts)  # - 1.
                twrp = twrp + xi * dsepv
                tbwrp = tbwrp + bxi * dsepv
            wrp[i] = 2.0 * twrp
            wrperr[i] = 2.0 * np.std(tbwrp) if nbts > 0 else 0.0
    return (wrp, wrperr)


# =============================================================================
def plotboot(cnt):
    """
    Just a test script to calculate and plot the w(rp) of each bootstrap sample,
    using Landy-Szalay estimator
    """
    npt, nrpt = 1.0 * cnt.npt, 1.0 * cnt.npt1
    nsepp, nsepv, dsepv, nbts = cnt.par.nsepp, cnt.par.nsepv, cnt.par.dsepv, cnt.par.nbts
    dd, rr, bdd = cnt.dd, cnt.rr, cnt.bdd
    nf = (nrpt / npt) * ((nrpt - 1) / (npt - 1))  # normalizing factor

    bxi = np.zeros(nbts)
    wrp = np.zeros(nsepp)
    wrperr = np.zeros(nsepp)
    wrpboot = np.zeros([nsepp, nbts])

    for i in range(nsepp):
        twrp = 0
        tbwrp = 0
        for j in range(nsepv):
            if rr[i, j] > 0.0:
                xi = nf * dd[i, j] / rr[i, j] - 1
                for k in range(nbts):
                    bxi[k] = nf * bdd[i, j, k] / rr[i, j] - 1
            else:
                xi = 0.0
                bxi = np.zeros(nbts)
            twrp = twrp + xi * dsepv
            tbwrp = tbwrp + bxi * dsepv
        wrp[i] = 2.0 * twrp  # 2.*
        wrpboot[i, :] = 2.0 * tbwrp
        wrperr[i] = 2.0 * np.std(tbwrp) if nbts > 0 else 0.0

    x = cnt.rpm
    for k in range(nbts):
        plotcf(x, wrpboot[:, k], yerr=0, color="g", marker=None, linewidth=0.2)
    plotcf(x, wrp, yerr=wrperr, color="r", marker=None, linewidth=1.5)

    return


# =============================================================================
def packpars(kind="pcf", **kwargs):
    """
    Pack a set input parameters into a **par** dictionary with attribute-style
    access (of class ``Munch``). Counting routines rely on this object to easily
    pass multiple parameters at once. For example:

    * par=packpars(kind='pcf') :
        Return default parameters, for a projected correlation run

    * par=packpars(kind='pcf',nsepp=12) :
        Return default parameters and set 12 bins in projected space

    The key-value pairs stored can be accessed with dot notation, e.g.
    ``par.nsepp``, ``par.omegam``, etc. and the quickest way to print all
    parameters nicely formatted is by typing :code:`par.qprint()`

    Any missing parameter in ``kwargs`` is set to its default value.
    After being created, you can add or modify values, e.g. ``par.nsepv=40``.
    See :ref:`in` for a detailed description of each parameter

    .. rubric :: Parameters

    kind : string. Default='pcf'
        * 'pcf' : Projected auto-correlation
        * 'pccf' : Projected cross-correlation
        * 'rcf' : Redshift space auto-correlation
        * 'rccf' : Resdhift space cross-correlation
        * 'acf' : Angular auto-correlation
        * 'accf' : Angular cross-correlation
        * 'rppiA' : Projected-radial space pair counts (auto-corr.)
        * 'rppiC' : Projected-radial space pair counts (cross-corr.)
        * 'sA' : Redshift space pair counts (auto-corr.)
        * 'sC' : Redshift space pair counts (cross-corr.)
        * 'thA' : Angular space pair counts (auto-corr.)
        * 'thC' : Angular space pair counts (cross-corr.)

    kwargs : [key]=[value]
        Optional keyword parameters

    .. rubric :: Returns

    par : Munch dictionary
        The **par** dictionary of input parameters

    .. rubric :: Examples

    .. code-block:: python

        par = packpars(kind='pcf', outfn='redgals', h0=70.)
        # All omitted parameters get their default value

        par.doboot = False
        # Any parameter can be set later with dot-style notation

        par.extra_irrelevant_info = 'yeah'
        # Even those you might want to invent on the fly

        par.qprint()
        # Prints a nicely formatted list of everything in par
    """
    p = Munch()
    p.kind = kind  # Type of par object
    p.autogrid = True  # Skip table autogrid
    p.dens = None  # Target density for SK autogrid
    p.outfn = "run001"  # Default path/base_name for output files
    p.pxorder = "natural"  # Pixel sorting method
    p.custRAbound = None  # Custom RA boundaries
    if kind == "thA" or kind == "thC":
        p.mxh1 = 20  # Nr of DEC bins for SK grid
        p.mxh2 = 20  # Nr of RA bins for SK grid
        p.doboot = False  # Do bootstrap counts
        p.nbts = 100  # Nr of boostrap samples
        p.bseed = 12345  # Seed for bootstrap weights
        p.wfib = False  # Apply fiber correction weights
        p.nsept = 36  # Nr of bins in theta [deg]
        p.septmin = 0.01  # Minimum theta [deg]
        p.dsept = 0.1  # Bin size in theta
        p.logsept = 1  # Do linear/log binning
        p.file = ""  # Filename of input data catalog
        p.description = ""  # Description of CF run
        p.cra = "ra"  # Colunm name of RA coordinate [deg]
        p.cdec = "dec"  # Colunm name of DEC coordinate [deg]
        p.cwei = "wei"  # Colunm name of source weight
    if kind == "thC":
        p.file1 = ""  # Filename of input cross catalog
        p.cra1 = "ra"  # Colunm name of RA coordinate in cross catalog [deg]
        p.cdec1 = "dec"  # Colunm name of DEC coordinate in cross catalog [deg]
        p.cwei1 = "wei"  # Colunm name of source weight in cross catalog
    if kind == "sA" or kind == "sC":
        p.h0 = 100.0  # Hubble constant [km/s/Mpc]
        p.omegam = 0.3  # Matter density
        p.omegal = 0.7  # Dark energy density
        p.mxh1 = 20  # Nr of DEC bins for SK grid
        p.mxh2 = 20  # Nr of RA bins for SK grid
        p.mxh3 = 20  # Nr of Z bins for SK grid
        p.doboot = False  # Do bootstrap counts
        p.nbts = 50  # Nr of boostrap samples
        p.bseed = 12345  # Seed for bootstrap weights
        p.wfib = False  # Apply fiber correction weights
        p.nseps = 36  # Nr of bins in redshift space s [Mpc/h]
        p.sepsmin = 0.01  # Minimum value of s [Mpc/h]
        p.dseps = 0.1  # Bin size in s
        p.logseps = 1  # Do linear/log binning
        p.calcdist = True  # Calculate comoving distances
        p.file = ""  # Filename of input data catalog
        p.description = ""  # Description of CF run
        p.cra = "ra"  # Colunm name of RA coordinate [deg]
        p.cdec = "dec"  # Colunm name of DEC coordinate [deg]
        p.cred = "z"  # Colunm name of REDSHIFT
        p.cwei = "wei"  # Colunm name of source weight
        p.cdcom = "dcom"  # Column name for comoving distance
    if kind == "sC":
        p.file1 = ""  # Filename of input cross catalog
        p.cra1 = "ra"  # Colunm name of RA coordinate in cross catalog [deg]
        p.cdec1 = "dec"  # Colunm name of DEC coordinate in cross catalog [deg]
        p.cred1 = "z"  # Colunm name of REDSHIFT in cross catalog
        p.cwei1 = "wei"  # Colunm name of source weight in cross catalog
        p.cdcom1 = "dcom"  # Column name for comoving distance
    if kind == "rppiA" or kind == "rppiC":
        p.h0 = 100.0  # Hubble constant [km/s/Mpc]
        p.omegam = 0.3  # Matter density
        p.omegal = 0.7  # Dark energy density
        p.mxh1 = 20  # Nr of DEC bins for SK grid
        p.mxh2 = 20  # Nr of RA bins for SK grid
        p.mxh3 = 20  # Nr of Z bins for SK grid
        p.doboot = False  # Do bootstrap counts
        p.nbts = 50  # Nr of boostrap samples
        p.bseed = 12345  # Seed for bootstrap weights
        p.wfib = False  # Apply fiber correction weights
        p.nsepp = 36  # Nr of bins in projected radius rp [Mpc/h]
        p.seppmin = 0.01  # Minimum value of rp [Mpc/h]
        p.dsepp = 0.1  # Bin size in rp
        p.logsepp = True  # Do linear/log binning in rp
        p.nsepv = 1  # Nr of bins in radial separation pi
        p.dsepv = 40.0  # Bin size in pi [Mpc/h]
        p.calcdist = True  # Calculate comoving distances
        p.file = ""  # Filename of input data catalog
        p.description = ""  # Description of CF run
        p.cra = "ra"  # Colunm name of RA coordinate [deg]
        p.cdec = "dec"  # Colunm name of DEC coordinate [deg]
        p.cred = "z"  # Colunm name of REDSHIFT
        p.cwei = "wei"  # Colunm name of source weight
        p.cdcom = "dcom"  # Column name for comoving distance
    if kind == "rppiC":
        p.file1 = ""  # Filename of input cross catalog
        p.cra1 = "ra"  # Colunm name of RA coordinate in cross catalog [deg]
        p.cdec1 = "dec"  # Colunm name of DEC coordinate in cross catalog [deg]
        p.cred1 = "z"  # Colunm name of REDSHIFT in cross catalog
        p.cwei1 = "wei"  # Colunm name of source weight in cross catalog
        p.cdcom1 = "dcom"  # Column name for comoving distance
    if kind == "pcf" or kind == "pccf":
        p.h0 = 100.0  # Hubble constant [km/s/Mpc]
        p.omegam = 0.3  # Matter density
        p.omegal = 0.7  # Dark energy density
        p.mxh1 = 30  # Nr of DEC bins for SK grid
        p.mxh2 = 180  # Nr of RA bins for SK grid
        p.mxh3 = 40  # Nr of Z bins for SK grid
        p.doboot = False  # Do bootstrap counts and errors
        p.nbts = 50  # Nr of boostrap samples
        p.bseed = 12345  # Seed for bootstrap weights
        p.wfib = False  # Apply fiber correction weights
        p.nsepp = 22  # Nr of bins in projected radius rp [Mpc/h]
        p.seppmin = 0.01  # Minimum value of rp [Mpc/h]
        p.dsepp = 0.15  # Bin size in rp
        p.logsepp = True  # Do linear/log binning in rp
        p.nsepv = 1  # Nr of bins in radial separation pi
        p.dsepv = 40.0  # Bin size in pi [Mpc/h]
        p.calcdist = True  # Calculate comoving distances
        p.file = ""  # Filename of input data catalog
        p.file1 = ""  # Filename of input random catalog
        p.description = ""  # Description of CF run
        p.estimator = "NAT"  # Estimator of CF
        p.cra = "ra"  # Colunm name of RA coordinate in data catalog [deg]
        p.cdec = "dec"  # Colunm name of DEC coordinate in data catalog [deg]
        p.cred = "z"  # Colunm name of REDSHIFT in data catalog
        p.cwei = "wei"  # Colunm name of source weight in data catalog
        p.cdcom = "dcom"  # Column name for comoving distance in data catalog
        p.cra1 = "ra"  # Colunm name of RA coordinate in random catalog [deg]
        p.cdec1 = "dec"  # Colunm name of DEC coordinate in random catalog [deg]
        p.cred1 = "z"  # Colunm name of REDSHIFT in random catalog
        p.cwei1 = "wei"  # Colunm name of source weight in random catalog
        p.cdcom1 = "dcom"  # Column name for comoving distance in random catalog
    if kind == "pccf":
        p.file2 = ""  # Filename of input cross catalog
        p.estimator = "DP"  # !Overwrite as only DP estimator is implemented
        p.cra2 = "ra"  # Colunm name of RA coordinate in cross catalog [deg]
        p.cdec2 = "dec"  # Colunm name of DEC coordinate in cross catalog [deg]
        p.cred2 = "z"  # Colunm name of redshift in cross catalog
        p.cwei2 = "wei"  # Colunm name of source weight in cross catalog
        p.cdcom2 = "dcom"  # Column name for comoving distance
    if kind == "rcf" or kind == "rccf":
        p.h0 = 100.0
        p.omegam = 0.3
        p.omegal = 0.7
        p.mxh1 = 20
        p.mxh2 = 20
        p.mxh3 = 20
        p.doboot = False
        p.nbts = 50
        p.bseed = 12345
        p.wfib = False
        p.nseps = 36
        p.sepsmin = 0.01
        p.dseps = 0.1
        p.logseps = 1
        p.calcdist = True  # Calculate comoving distances
        p.file = ""
        p.file1 = ""
        p.description = ""
        p.estimator = "NAT"
        p.cra = "ra"
        p.cdec = "dec"
        p.cred = "z"
        p.cwei = "wei"
        p.cdcom = "dcom"  # Column name for comoving distance
        p.cra1 = "ra"
        p.cdec1 = "dec"
        p.cred1 = "z"
        p.cwei1 = "wei"
        p.cdcom1 = "dcom"  # Column name for comoving distance
    if kind == "rccf":
        p.estimator = "DP"  # only DP estimador is implemented
        p.cra2 = "ra"
        p.cdec2 = "dec"
        p.cred2 = "z"
        p.cwei2 = "wei"
        p.cdcom2 = "dcom"  # Column name for comoving distance
        p.file2 = ""
    if kind == "acf" or kind == "accf":
        p.mxh1 = 20
        p.mxh2 = 20
        p.doboot = False
        p.nbts = 50
        p.bseed = 12345
        p.wfib = False
        p.nsept = 36
        p.septmin = 0.01
        p.dsept = 0.1
        p.logsept = 1
        p.file = ""
        p.file1 = ""
        p.description = ""
        p.estimator = "NAT"
        p.cra = "ra"
        p.cdec = "dec"
        p.cwei = "wei"
        p.cra1 = "ra"
        p.cdec1 = "dec"
        p.cwei1 = "wei"
    if kind == "accf":
        p.estimator = "DP"  # only DP estimator is implemented
        p.cra2 = "ra"
        p.cdec2 = "dec"
        p.cred2 = "z"
        p.cwei2 = "wei"
        p.file2 = ""
    for key, value in kwargs.items():
        p.update({key: value})
    if "outfn" not in kwargs:
        print("Missing outfn: default to", p.outfn)

    return p


# =============================================================================
def radec2xyz(ra, dec, r=0.5):
    """
    Converts a (ra,dec) coordinates to rectangular (x,y,z) coordinates in a
    sphere of radius r.

    For ``r=0.5``, this allows to speed up subsecuent haversine distance
    calculations between two points, by simply computing
    :math:`dhav^2 = (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2`

    .. rubric :: Parameters

    ra,dec : float arrays
        Right ascention an declination coordinates
    r : float. Default=0.5
        Sphere radius

    .. rubric :: Returns

    (x,y,z) : tuple
        Rectangular coordinates
    """
    r = 0.5
    return (r * np.cos(dec) * np.sin(ra), r * np.cos(dec) * np.cos(ra), r * np.sin(dec))


# =============================================================================
def bound3d(decs, dcs):
    """
    Return the (maximal) survey boundaries in (RA,DEC,DCOM) for multiple samples.
    Note in the RA direction, the limits are **always** set as ramin=0. and
    ramax=360. to make the pair counting valid for surveys that cross the origin
    of coordinates

    .. rubric :: Parameters

    decs : float array or list of arrays
        DEC coordinates of one or more samples [deg]
    dcs : float array or list of arrays
        Comoving distances of one or more samples [Mpc/h]

    .. rubric :: Returns

    bound : tuple
        The limits (ramin,ramax,decmin,decmax,dcmin,dcmax) that enclose all
        input samples
    """
    delta = 0.001  # small delta to avoid issues with souces exactly at edges
    if type(decs) != list:
        decs = [decs]
        dcs = [dcs]
    ramin = 0.0
    ramax = 360.0
    decmin = max(min([min(d) for d in decs]) - delta, -90.0)
    decmax = min(max([max(d) for d in decs]) + delta, 90.0)
    dcmin = max(min([min(dc) for dc in dcs]) - delta, 0.0)
    dcmax = max([max(dc) for dc in dcs]) + delta
    return (ramin, ramax, decmin, decmax, dcmin, dcmax)


# =============================================================================
def bound2d(decs):
    """
    Return the (maximal) survey boundaries in (RA,DEC) for multiple samples.
    Note in the RA direction, the limits are **always** set as ramin=0. and
    ramax=360. to make the pair counting valid for surveys that cross the origin
    of coordinates

    .. rubric :: Parameters

    decs : float array or list of arrays
        DEC coordinates of one or more samples [deg]

    .. rubric :: Returns

    bound : tuple
        The limits (ramin,ramax,decmin,decmax) that enclose all input samples
    """
    delta = 0.001  # small delta to avoid issues with souces exactly at edges
    if type(decs) != list:
        decs = [decs]
    ramin = 0.0
    ramax = 360.0
    decmin = max(min([min(d) for d in decs]) - delta, -90.0)
    decmax = min(max([max(d) for d in decs]) + delta, 90.0)
    return (ramin, ramax, decmin, decmax)


# =============================================================================
def cross0guess(ra):
    """
    Guess if a set of RA coordinates cross the RA=0 division (by finding one
    source 1deg to the left and another 1deg to the right of RA=0, at least)

    .. rubric :: Parameters

    ra : array
        Right ascention coordiantes

    .. rubric :: Returns

    res : bool
        True if the sample seems to cross the RA=0 boundary
    """
    res = (np.where(ra > 359.0)[0].size) > 0 and (np.where(ra < 1.0)[0].size) > 0
    return res


# =============================================================================
def logtimming(log, cntid, t):
    """
    Write a message to a log instance, showing the compute time for the counts
    identified with a certain ID string

    .. rubric :: Parameters

    log : logging object
        Log object
    cntid : string
        String to ID a certain type of counts, e.g. 'DD', 'RR', DR', etc.
    t : float
        The time elapsed [s]
    """
    log.info(f"Done with {cntid} calculations. Looptime (s) : {t:0.3f}")
    sys.stdout.flush()


# =============================================================================
def setlogs(par, runspyder=True):
    """
    Set up the log machinery used by the main counting routines by creating the
    logger object, adding the required handlers and cleaning previous logs if
    present

    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters dictionary for Gundam routines
    runspyder : bool. Default=True
        If ``runspyder=True``, add an extra handler for stdout when running
        under a Spyder console

    .. rubric :: Returns

    log : logging object
        The log object
    logfile : string
        The complete path of the .log file
    """
    log = getLogger("cf")  # fix cf name for now
    logfile = os.path.abspath(par.outfn + ".log")
    logfileH = FileHandler(logfile, "w")
    log.addHandler(logfileH)
    if runspyder:
        log.addHandler(StreamHandler(sys.stdout))  # add only under Spyder
    log.setLevel(INFO)

    # Delete the Fortran log file (to avoid being appended in succesive runs)
    logfort = par.outfn + ".fortran.log"
    if os.path.isfile(logfort):
        os.remove(logfort)

    return (log, logfile)


# =============================================================================
def closelog(log, runspyder=True):
    """
    Close the log machinery used by the main counting routines by removing the
    handlers

    .. rubric :: Parameters

    log : logging object
        The log object
    runspyder : bool. Default=True
        If ``runspyder=True``, remove the extra handler for stdout added when
        running under a Spyder console
    """
    # this only works if handler[0] is a file handler
    if runspyder:
        log.removeHandler(log.handlers[1])  # remove only under Spyder
    log.handlers[0].close()
    log.removeHandler(log.handlers[0])


# =============================================================================
def check_kind(par, kind):
    """
    Check that ``par.kind = kind``. Useful to test, for example, if you are
    passing the right **par** object before really counting pairs

    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters dictionary for Gundam routines
    kind : string
        The `kind` to check against, e.g. 'pcf', 'accf', 'rppiA', etc.
    """
    msg = 'This function needs a parameter bject with par.kind = "' + kind + '"'
    if par.kind != kind:
        raise Exception(msg)


# =============================================================================
def logcallinfo(log, par, npts=[]):
    """
    Write to log some useful runtime parameters of the main counting routines

    .. rubric :: Parameters

    log : logging object
        The log object
    par : Munch dictionary
        Input parameters dictionary for Gundam routines
    npts : list of 1, 2 or 3 integers
        The nr of objects in the input data table, random table, and cross sample
        table
    """
    log.info("Calling Information  ================================================")
    if len(npts) > 0:
        msg = "table        = " + par.file
        log.info(msg)
        msg = "table_nrows  = " + str(npts[0])
        log.info(msg)
    if len(npts) > 1:
        msg = "table1       = " + par.file1
        log.info(msg)
        msg = "table1_nrows = " + str(npts[1])
        log.info(msg)
    if len(npts) > 2:
        msg = "table2       = " + par.file2
        log.info(msg)
        msg = "table2_nrows = " + str(npts[2])
        log.info(msg)
    msg = "nthreads     = " + str(par.nthreads)
    log.info(msg)
    msg = "write        = " + str(par.write)
    log.info(msg)
    msg = "plot         = " + str(par.plot)
    log.info(msg)
    log.info("=====================================================================")


# =============================================================================
def initialize(kind, par, nthreads=None, write=None, plot=None):
    """
    Perform some initialization tasks common to all correlation function runs,
    such as reading **par** from disk file if needed, check the parameter `kind`,
    find out if running under Spyder, initialize the logs, etc.

    .. rubric :: Parameters

    kind : string
        The `kind` to check against, e.g. 'pcf', 'accf', 'rppiA', etc.
    par : Munch dictionary
        Input parameters dictionary for Gundam routines
    nthreads : integer.
        Number of threads to use. Passed just to store it under **par**
    write : bool. Default=True
        Flag to generate output count files. Passed just to store it under **par**
    plot : bool. Default=True
        Flag to generate plot. Passed just to store it under **par**

    .. rubric :: Returns

    par : Munch dictionary
        Input parameters dictionary for Gundam routines
    log : logging object
        The log object
    logfile, logfilefort : string
        The complete path of the `.log` file and the `.fotran.log` file
    runspyder : bool. Default=True
        Detect if running under a Spyder console
    """
    # Load parameter object from file
    if type(par) is str:
        par = readpars(par)

    # Check we have the right kind parameters
    check_kind(par, kind)

    # Find out if running under Spyder via runfile() or directly in a terminal
    # to avoid double output while running in terminal
    runspyder = False if sys.stdin.isatty() else True

    # Set the log files  ------------------------------------------------------
    log, logfile = setlogs(par, runspyder=runspyder)
    logfilefort = par.outfn + ".fortran.log"

    # Log the start time and output file name  --------------------------------
    t0 = time.time()
    log.info(time.strftime("START!   %x  %H:%M:%S"))
    log.info("Name of output set to : " + par.outfn)

    # Add runtime call information to parameter object  -----------------------
    par.nthreads = nthreads
    par.write = write
    par.plot = plot

    return (par, log, logfile, logfilefort, runspyder, t0)


# =============================================================================
def finalize(log, logf, logff, Ltime, t0, counts):
    """
    Perform some finalization tasks common to all correlation runs, by logging
    loop/total times and adding the contents of log files to the **counts**
    output dictionary

    .. rubric :: Parameters

    log : logging object
        The log object
    logfile : string
        The complete path of the `.log` file
    Ltime, Ttime : float
        Loop time and total compute time of a correlation run
    counts : Munch dictionary
        Output dictionary containing all counts and correlations
    """
    t1 = time.time()
    log.info(time.strftime("DONE!   %x  %H:%M:%S"))
    log.info("Loop | Total time (s) ".ljust(27) + f" : {Ltime:0.2f} | {t1-t0:0.2f}")
    addlog(logf, "log", counts)
    addlog(logff, "logfortran", counts)


# =============================================================================
def tidy_counts(tt, par):
    """
    Do some tiding of count/boostrap counts arrays returned by external Fortran
    routines. Basically: (1) Separate dd/bdd, if present, and (2) transpose arrays
    from Fortran optimized indexing back to a more natural row-oriented indexing

    .. rubric :: Parameters

    tt : list of ndarray
        Ouput counts as returned by Fortran counting routines
    par : Munch dictionary
        Input parameters dictionary

    .. rubric :: Returns

    (dd,bdd) : tuple of arrays
        Count and boostrap count arrays. If ``doboot=False``, ``bdd`` is set to a
        zero-valued array with its boostrap samples dimension zeroed
    """
    if par.doboot:
        dd, bdd = tt
    else:
        dd = tt
        if par.kind in ["pcf", "pccf", "rppiA", "rppiC"]:
            bdd = np.zeros([0, par.nsepv, par.nsepp])
        if par.kind in ["acf", "accf", "thA", "thC"]:
            bdd = np.zeros([0, par.nsept])
        if par.kind in ["rcf", "rccf", "sA", "sC"]:
            bdd = np.zeros([0, par.nseps])

    # Transpose resuls back to a more natural order in python, as Fortran
    # routines return arrays with dimensions exchanged to optimize cache utilization
    if par.kind in ["pcf", "pccf", "rppiA", "rppiC"]:
        dd = dd.transpose([1, 0])  # orig=[1,0]    invert=[0,1]
        bdd = bdd.transpose([1, 2, 0])  # orig=[2,1,0]  invert=[0,1,2]
    if par.kind in ["acf", "accf", "rcf", "rccf", "thA", "sA", "thC", "sC"]:
        bdd = bdd.transpose([1, 0])
    return (dd, bdd)


# =============================================================================
def addlog(file, key, m):
    """
    Read a text file and dump it into a key of a given dictionary. Useful to
    add entire logs from disk into Munch objects

    .. rubric :: Parameters

    file : string
        Complete path of file, usually any log file
    key : string
        Name of the key to be created
    m : Munch dictionary
        The dictionary where the key ``m.key`` will be created
    """
    try:
        with open(file, "r") as f:  # save log contents
            contents = f.read()
            m[key] = contents
    except:
        pass


# =============================================================================
def pixsort(tab, colnms, par):
    """
    Arrange an astropy table with (ra,dec) or (ra,dec,z) data into a grid of
    SK pixels and sort the rows of the table according to the pixel index, ordered
    by a given methodology. This way data in a given pixel sits closer in memory,
    largely increasing the efficiency of the cache.

    Several (experimental) orders such as Morton and Hilbert are available

    .. rubric :: Parameters

    tab: astropy table
        Data table with ra/dec or ra/dec/z points
    colnms: list of strings
        Colum names for coordinates and redshift, e.g. ['ra','dec','z'],
        ['RAJ2000','DEJ2000','redshift']
    par: Munch dictionary
        Input parameters dictionary. Must contain ``mxh1``, ``mxh2``, ``mxh3``
        (when appropiate), the survey boundaries ``sbound`` and the desired method
        for ordering ``pxorder``. For samples that cross the RA=0 limit, it can
        is useful to specify a custom RA boundary. See :ref:`Custom RA boundaries <custRAbound>`

    .. rubric :: Returns

    sidx : array
        Sort index. Use ``tab=tab[sidx]`` to actually sort the data
    """
    if par.kind in ["pcf", "pccf", "rcf", "rccf", "rppiA", "rppiC", "sA", "sC"]:
        cra, cdec, cred = colnms
        ra = tab[cra].data
        zmin, zmax = min(tab[cred].data), max(tab[cred].data)
        ramin, ramax, decmin, decmax, dcmin, dcmax = par.sbound
        if par.custRAbound is not None:
            r1 = 360.0 - par.custRAbound[0]
            ra += r1  # in place addition to avoid copy // ra = ra + r1
            ra[ra >= 360.0] = ra[ra >= 360.0] - 360.0
            ramin = 0.0
            ramax = par.custRAbound[1] + r1
            binsra = np.linspace(ramin, ramax, num=par.mxh2)
        else:
            binsra = np.linspace(ramin, ramax, num=par.mxh2)
        pxra = np.digitize(tab[cra].data, bins=binsra)  # opt. add as column tab['pxra']
        binsdec = np.linspace(decmin, decmax, num=par.mxh1)
        pxdec = np.digitize(tab[cdec].data, bins=binsdec)  # opt. add as column tab['pxdec']
        binsred = np.linspace(zmin, zmax, num=np.int_(par.mxh3))  # num=par.mxh3*10
        pxred = np.digitize(tab[cred].data, bins=binsred)  # opt. add as column tab['pxred']

        # CHOOSE THE ORDERING METHOD  ================================
        # My own "pixel sorting" where data is ordered according to its
        # corresponding pixel index in DEC, RA and Z, succesively.
        if par.pxorder == "natural":
            # lexsort() is ~2x faster than argsort(), e.g. sidx=tab.argsort(['pxdec','pxra','pxred'])
            # Note sorting keys must be given in reverse order, i.e. from last to first
            sidx = np.lexsort((pxra, pxdec, pxred))
            # sidx = np.lexsort((pxred, pxra, pxdec))

        # Morton order in RA-DEC, using the RA-DEC FP values multiplied by a
        # large(?) integer
        elif par.pxorder == "morton_dir":
            import mortcode as mc

            lon = tab["RA"].data
            lat = tab["DEC"].data
            if min(lat) < 0.0:
                lat = lat - min(lat)
            mort = list(map(mc.get_latlong_morton, lat, lon))
            sidx = np.argsort(mort)
            tab["mort"] = mort

        # Morton order in RA-DEC, using directly the pixel values to sort
        elif par.pxorder == "morton2":
            import mortcode as mc

            mort = list(map(mc.get_morton, pxra, pxdec))
            sidx = np.argsort(mort)
            tab["mort"] = mort

        # Morton order in RA-DEC-Z, using directly the pixel values to sort
        elif par.pxorder == "morton3":
            import pymorton as pm

            mort = list(map(pm.interleave3, pxra, pxdec, pxred))
            sidx = np.argsort(mort)
            tab["mort"] = mort

        # Hilbert order in RA-DEC, using directly the pixel values to sort
        elif par.pxorder == "hilbert2":
            from hilbert_curve import xy2d

            order = 8
            hilb = list(map(xy2d, [order] * len(pxra), pxra, pxdec))
            sidx = np.argsort(hilb)
            tab["hilb"] = hilb

        # Hilbert order in RA-DEC-Z, using directly the pixel values to sort
        elif par.pxorder == "hilbert3":
            from ndhilbert import Hilbert

            h = Hilbert(3)
            hilb = list(map(h.decode, list(zip(pxra, pxdec, pxred, strict=False))))
            sidx = np.argsort(hilb)
            tab["hilb"] = hilb

        else:
            raise NameError("Order " + par.pxorder + " not implemented")

        # Optionally add pixel, morton, hilbert numbers to data table. Nice for plotting
        # tab['pxdec'] = pxdec
        # tab['pxra']  = pxra
        # tab['pxred'] = pxred
    else:
        cra, cdec = colnms
        ra = tab[cra].data
        ramin, ramax, decmin, decmax = par.sbound
        if par.custRAbound is not None:
            r1 = 360.0 - par.custRAbound[0]
            ra += r1  # in place addition to avoid copy // ra = ra + r1
            ra[ra >= 360.0] = ra[ra >= 360.0] - 360.0
            ramin = 0.0
            ramax = par.custRAbound[1] + r1
            binsra = np.linspace(ramin, ramax, num=par.mxh2)
        else:
            binsra = np.linspace(ramin, ramax, num=par.mxh2)
        pxra = np.digitize(tab[cra].data, bins=binsra)  # opt. add as column tab['pxra']
        binsdec = np.linspace(decmin, decmax, par.mxh1)
        pxdec = np.digitize(tab[cdec].data, bins=binsdec)  # opt. add as column tab['pxdec']

        # CHOOSE THE ORDERING METHOD  ================================
        # My own "pixel sorting" where data is ordered according to its
        # corresponding pixel index in DEC, RA and Z, succesively.
        if par.pxorder == "natural":
            # lexsort() is ~2x faster than argsort(), e.g. sidx=tab.argsort(['pxdec','pxra','pxred'])
            # Note sorting keys must be given in reverse order, i.e. from last to first
            sidx = np.lexsort((pxra, pxdec))

        # Morton order in RA-DEC, using directly the pixel values to sort
        elif par.pxorder == "morton2":
            import mortcode as mc

            mort = list(map(mc.get_morton, pxra, pxdec))
            sidx = np.argsort(mort)
            tab["mort"] = mort

        elif par.pxorder in [None, "none", "None"]:
            sidx = np.arange(0, len(tab))

        else:
            raise NameError("Order " + par.pxorder + " not implemented")

        # Optionally add pixel, morton, hilbert numbers to data table. Nice for plotting
        # tab['pxdec'] = pxdec
        # tab['pxra']  = pxra

    return sidx


# =============================================================================
def bestSKgrid2d(par, npts, ras, dens=None):
    """
    Try to find the optimum size (mxh1,mxh2) of a 2D skip grid, for an arbitrary
    sample of (ra,dec) points.

    This is far from trivial, so for the moment, this routine works as follows:

    #. Choose an "optimum" target cell density (around 22 ppcell in many tests)

    #. Find the best ``mxh1`` from an empirical best fit relation of npts vs mxh1

    #. Adjust ``mxh2`` to reach the target cell density


    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters dictionary
    npts : integer or list of integers
        Nr. of objects in the sample(s). For cross-correlations it should
        be a list of the form [sample_D_size, sample_R,size] and
        the effective ``npts`` adopted is their sum
    ras : float array or list of arrays
        Right ascention of objects in the samples(s). For cross-correlations,
        the two samples are concatenated
    dens : float. Default=None
        Target cell density. Note that `dens=None` implicates a default value
        of 22 if npts>50000 and 16 otherwise.

    .. rubric :: Returns

    mxh1, mxh2 : integer
        Nr. of DEC and RA cells of the SK table, respectively
    dens : float
        The effective cell density adopted
    """
    if par.kind in ["thC"]:
        ras = np.concatenate(ras)  # Flatten the list of two arrays if needed
        npts = sum(npts)  # In DR case, sum the size of both samples

    # Set default density based on the best nr.part.per.cell as function of
    # sample size, as found in timming grids.
    if dens is None:
        if par.kind in ["thA", "thC"]:  # potentially tune for thC case
            dens = 22.0 if npts > 50000 else 16.0  # best times have sort of constant density
    # print('SK cell target density : ', dens)

    # Find real limits in RA among all input samples
    ramin, ramax = min(ras), max(ras)
    if par.custRAbound is not None:
        samplewidth = (360 - par.custRAbound[0]) + par.custRAbound[1]
    else:
        samplewidth = ramax - ramin

    h1h2 = npts / dens  # Combined h1h2

    # Choose mxh1  ------------------------------
    if par.kind in ["thA", "thC"]:  # potentially tune for thC case
        # Based on the fitting of a+b*sqrt(n) to n vs best_mhx1 for thA case
        h1 = max(np.int_(np.rint(10.75 + 0.075 * np.sqrt(npts))), 1)

    # Choose mxh2  ------------------------------
    # Set to reach the target density (last factor is to extend to 360 range)
    h2 = np.int_(np.rint(h1h2 / h1) * (360.0 / samplewidth))

    # Implement some (min,max) safeguards in h1,h2 ?
    return [h1, h2, dens]


# =============================================================================
def bestSKgrid3d(par, npts, ras, dens=None):
    """
    Try to find the optimum size (mxh1,mxh2,mxh3) of a 3D skip grid, for an
    arbitrary sample of (ra,dec,dcom) points.

    This is far from trivial, so for the moment, this routine works as follows:

    #. Choose an "optimum" target cell density according to the type of correlation

    #. Choose the best ``mxh3`` as ``mxh3=int((dcmax-dcmin)/rvmax)``

    #. Find the best ``mxh1`` from an empirical best fit relation of npts vs mxh1

    #. Adjust ``mxh2`` to reach the target cell density


    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters dictionary
    npts : integer or list of integers
        Nr. of objects in the sample(s). For cross-correlations it should
        be a list of the form [sample_D_size, sample_R,size] and
        the effective ``npts`` adopted is their sum
    ras : float array or list of arrays
        Right ascention of objects in the samples(s). For cross-correlations,
        the two samples are joined together
    dens : float. Default=None
        Target cell density. Note that `dens=None` implicates different default
        values. See function code.


    .. rubric :: Returns

    mxh1, mxh2, mxh3 : integer
        Nr. of DEC, RA, DCOM cells of the SK table, respectively
    dens : float
        The effective cell density adopted
    """
    if par.kind in ["rppiC", "sC"]:
        ras = np.concatenate(ras)  # Flatten the list of two arrays if needed
        npts = sum(npts)  # In DR case, sum the size of both samples

    # Set default density based on the best nr.part.per.cell as function of
    # sample size, as found in timming grids.
    if dens is None:
        if par.kind in ["rppiA", "rppiC"]:
            dens = 18.0 if npts > 100000 else 8.0  # best times have sort of constant density
        elif par.kind in ["sA", "sC"]:
            dens = 28.0 if npts > 100000 else 12.0  # best times have sort of constant density

    # print('SK cell target density : ', dens)

    # Find real limits among all input samples
    dummramin, dummramax, decmin, decmax, dcmin, dcmax = par.sbound
    ramin, ramax = min(ras), max(ras)
    if par.custRAbound is not None:
        samplewidth = (360 - par.custRAbound[0]) + par.custRAbound[1]
    else:
        samplewidth = ramax - ramin

    # Choose mxh3  ------------------------------
    if par.kind in ["rppiA", "rppiC"]:
        radmax = makebins(par.nsepv, 0.0, par.dsepv, 0)[0][-1]
    elif par.kind in ["sA", "sC"]:
        radmax = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)[0][-1]
    h3 = np.int_((dcmax - dcmin) / radmax)

    h1h2 = npts / (dens * h3)  # Combined h1h2

    # Choose mxh1  ------------------------------
    if par.kind in ["rppiA", "rppiC"]:
        # Based on the fitting of a+b*sqrt(n) to n vs best_mhx1 of rrpiA
        h1 = max(np.int_(np.rint(2.92 + 0.05 * np.sqrt(npts))), 1)
    elif par.kind in ["sA", "sC"]:
        # Based on the fitting of a+b*sqrt(n) to n vs best_mhx1 if sA
        h1 = max(np.int_(np.rint(4.03 + 0.03 * np.sqrt(npts))), 1)

    # Choose mxh2  ------------------------------
    # Set to reach the target density (last factor is to extend to 360 range)
    h2 = np.int_(np.rint(h1h2 / h1) * (360.0 / samplewidth))

    # Implement some (min,max) safeguards in h1,h2,h3 ?
    return [h1, h2, h3, dens]


# =============================================================================
def pairs_auto(par, wunit, logff, tab, x, y, z, sk, ll, dc=None):
    """
    Wrapper for calling Fortran counting routines (for auto-pairs)

    This function isolates the call of external Fortran routines, choosing
    the correct one based on geometry, and choosing the fastest depending whether,
    weights, bootstrap errors, etc. are requested or not.


    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters

    wunit : bool
       * True : all weights are equal to 1
       * False : at least one weight is not 1

    logff : string
        File for logging the output of Fortran routines

    tab : astropy table
        Table with data particles

    x,y,z : arrays
        Rectangular coordinates of data particles

    sk,ll : arrays
        Skip table (SK) and linked list table (see :func:`gundam.skll2d` or :func:`gundam.skll3d`)

    dc : array [optional]
        Array of comoving distances. Not needed for angular counts


    .. rubric :: Returns

    tt : list of ndarray
        Ouput counts as returned by Fortran counting routines
    """
    nt = par.nthreads
    npt = len(tab)

    if par.kind in ("rppiA"):  ############  Proyected Space  ########### 'pcf'
        sepp = makebins(par.nsepp, par.seppmin, par.dsepp, par.logsepp)[0]
        sepv = makebins(par.nsepv, 0.0, par.dsepv, False)[0]
        if par.doboot:  # Counting Data + Boostrap
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.rppi_Ab(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.rppi_Ab_wg(*args)  # weighted counting
        else:  # Counting Data only
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.rppi_A(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.rppi_A_wg(*args)  # weighted counting

    if par.kind in ("sA"):  ###########  Redsfhift Space   ########### 'rcf'
        seps = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)[0]
        if par.doboot:  # Counting Data + Boostrap
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.s_Ab(*args)  # fast unweighted counting
            else:
                # args = [nt,npt,tab[par.cdec].data,dc,tab[par.cwei],x,y,z,par.nsepp,sepp,par.nsepv,sepv,par.sbound,par.mxh1,par.mxh2,par.mxh3,par.nbts,par.bseed,par.wfib,par.cntid,logff,sk,ll]
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.s_Ab_wg(*args)  # weighted counting
        else:  # Counting Data only
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.s_A(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.s_A_wg(*args)  # weighted counting

    if par.kind in ("thA"):  #############  Angular Space  ############## 'acf'
        sept = makebins(par.nsept, par.septmin, par.dsept, par.logsept)[0]
        if par.doboot:  # Counting Data + Boostrap
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    x,
                    y,
                    z,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.nbts,
                    par.bseed,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.th_Ab(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.nbts,
                    par.bseed,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.th_Ab_wg(*args)  # weighted counting
        else:  # Counting Data only
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    x,
                    y,
                    z,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.th_A(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cdec].data,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk,
                    ll,
                ]
                tt = cff.mod.th_A_wg(*args)  # weighted counting

    return tt


# =============================================================================
def pairs_cross(par, wunit, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1, dc=None, dc1=None):
    """
    Wrapper for calling Fortran counting routines (for cross-pairs among two
    tables called D(data) and Random(R))

    This function isolates the call of external Fortran routines, choosing
    the correct one based on geometry, and choosing the fastest depending whether,
    weights, bootstrap errors, etc. are requested or not.


    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters

    wunit : bool
       * True : all weights are equal to 1
       * False : at least one weight is not 1

    logff : string
        File for logging the output of Fortran routines

    tab,tab1 : astropy tables
        Table with data particles (D) and random particles (R), respectively

    x,y,z,x1,y1,z1 : arrays
        Rectangular coordinates of particles in D table and R rable, respectively

    tab1 : astropy table
        Table with data particles (D)

    sk1,ll1 : arrays
        Skip table (SK) and linked list table (see :func:`gundam.skll2d` or :func:`gundam.skll3d`)

    dc,dc1 : arrays [optional]
        Array of comoving distances of particles in D and R tables. Not needed for angular counts


    .. rubric :: Returns

    tt : list of ndarray
        Ouput counts as returned by Fortran counting routines
    """
    nt = par.nthreads
    npt = len(tab)
    npt1 = len(tab1)

    if par.kind in ("rppiC"):  ###########  Proyected Space  ########### 'pccf'
        sepp = makebins(par.nsepp, par.seppmin, par.dsepp, par.logsepp)[0]
        sepv = makebins(par.nsepv, 0.0, par.dsepv, False)[0]
        if par.doboot:  # Counting Data + Boostrap
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    x1,
                    y1,
                    z1,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.rppi_Cb(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei].data,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    tab1[par.cwei1].data,
                    x1,
                    y1,
                    z1,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.rppi_Cb_wg(*args)  # weighted counting
        else:  # Counting Data only
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    x1,
                    y1,
                    z1,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.rppi_C(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei].data,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    tab1[par.cwei1].data,
                    x1,
                    y1,
                    z1,
                    par.nsepp,
                    sepp,
                    par.nsepv,
                    sepv,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.rppi_C_wg(*args)  # weighted counting

    if par.kind in ("sC"):  ############  Proyected Space  ############# 'rccf'
        seps = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)[0]
        if par.doboot:  # Counting Data + Boostrap
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    x1,
                    y1,
                    z1,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.s_Cb(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei].data,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    tab1[par.cwei1].data,
                    x1,
                    y1,
                    z1,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.nbts,
                    par.bseed,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.s_Cb_wg(*args)  # weighted counting
        else:  # Counting Data only
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    x1,
                    y1,
                    z1,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.s_C(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    dc,
                    tab[par.cwei].data,
                    x,
                    y,
                    z,
                    npt1,
                    dc1,
                    tab1[par.cwei1].data,
                    x1,
                    y1,
                    z1,
                    par.nseps,
                    seps,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.mxh3,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.s_C_wg(*args)  # weighted counting

    if par.kind in ("thC"):  ############  Angular Space  ############## 'accf'
        sept = makebins(par.nsept, par.septmin, par.dsept, par.logsept)[0]
        if par.doboot:  # Counting Data + Boostrap
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    x,
                    y,
                    z,
                    npt1,
                    x1,
                    y1,
                    z1,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.nbts,
                    par.bseed,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.th_Cb(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    tab[par.cwei].data,
                    x,
                    y,
                    z,
                    npt1,
                    tab1[par.cwei1].data,
                    x1,
                    y1,
                    z1,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.nbts,
                    par.bseed,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.th_Cb_wg(*args)  # weighted counting
        else:  # Counting Data only
            if par.wfib == False and wunit:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    x,
                    y,
                    z,
                    npt1,
                    x1,
                    y1,
                    z1,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.th_C(*args)  # fast unweighted counting
            else:
                args = [
                    nt,
                    npt,
                    tab[par.cra].data,
                    tab[par.cdec].data,
                    tab[par.cwei],
                    x,
                    y,
                    z,
                    npt1,
                    tab1[par.cwei1],
                    x1,
                    y1,
                    z1,
                    par.nsept,
                    sept,
                    par.sbound,
                    par.mxh1,
                    par.mxh2,
                    par.wfib,
                    par.cntid,
                    logff,
                    sk1,
                    ll1,
                ]
                tt = cff.mod.th_C_wg(*args)  # weighted counting

    return tt


# =============================================================================
def rppi_A(tab, par, nthreads=-1, write=True, plot=False, **kwargs):
    r"""
    Given an astropy data table, count pairs in the projected-radial space
    (:math:`r_p` vs :math:`\pi` space).

    All input parameters that control binning, cosmology, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles

    par : Munch dictionary
        Input parameters. See :ref:`indic-rppiA` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=True
       * True : generate a plot of counts (integrated radially) vs proj. radius
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dd, counts.bdd, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicrppiA`
        for a detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                     # read data
        par = gun.packpars(kind='rppiA',outfn='redgalpairs')  # generate default parameters
        cnt = gun.rppiA(gals, par, write=True, plot=True)     # get pair counts and plot
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF funcs
    (par, log, logf, logff, runspyder, t0) = initialize(
        "rppiA", par, nthreads=nthreads, write=write, plot=plot
    )

    # Find number of particles in input table and set nr of threads  ----------
    npt = len(tab)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt])

    # Create bins in p(projected) and v(line-of-sight) space ------------------
    sepp, seppout = makebins(par.nsepp, par.seppmin, par.dsepp, par.logsepp)
    sepv, sepvout = makebins(par.nsepv, 0.0, par.dsepv, False)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data

    # Define the boundaries of the sample  ------------------------------------
    par.sbound = bound3d(tab[cdec].data, dc)
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0".ljust(lj) + " : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # If requested, try to find the best SK grid size  ------------------------
    if par.autogrid:
        log.info("SK Autogrid".ljust(lj) + " : ON")
        par.mxh1, par.mxh2, par.mxh3, tdens = bestSKgrid3d(par, npt, tab[cra].data, dens=par.dens)
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("SK Autogrid".ljust(lj) + " : OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par.mxh1, par.mxh2, par.mxh3]))

    # Sort data table/s according to some order  ------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par)
    tab, dc = tab[sidx], dc[sidx]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  ------------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt, tab[cra], tab[cdec], dc, par.sbound, sepv, par.nsepv
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------------------------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1.0).all()

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    par.cntid = "DD"
    log.info("====  Counting " + par.cntid + " pairs in " + str(par.mxh1) + " DEC bands  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_auto(par, wunit, logff, tab, x, y, z, sk, ll, dc=dc)
    tend = time.time()
    logtimming(log, par.cntid, tend - tstart)
    tacc = tend - tstart
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts and integrate along pi direction  ---------------------
    dd, bdd = tidy_counts(tt1, par)
    intpi = dd.sum(axis=1)
    intpib = bdd.sum(axis=1) if par.doboot else None

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("Counts plot")
            plotcf(seppout[1], intpi, 0, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutputC(par, npts=[npt], binslmr=seppout, dd=dd, bootc=bdd, intpi=intpi, intpib=intpib)

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_rppicounts(*seppout, dd, par, cntid="rppi")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def rppi_C(tab, tab1, par, nthreads=-1, write=True, plot=False, **kwargs):
    r"""
    Given two astropy data tables, cross-count pairs in the projected-radial space
    (:math:`r_p` vs :math:`\pi` space).

    All input parameters that control binning, cosmology, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with particles (sample D)

    tab1 : astropy table
        Table with particles (sample R)

    par : Munch dictionary
        Input parameters. See :ref:`indic-rppiC` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=True
       * True : generate a plot of counts (integrated radially) vs proj. radius
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dr, counts.bdr, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicrppiC`
        for a detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        qsos = Table.read('qso.fits')                          # read data
        gals = Table.read('redgal.fits')                       # read data
        par = gun.packpars(kind='rppiC',outfn='qso_rg_pairs')  # generate default parameters
        cnt = gun.rppiC(qsos, gals, par, write=True)           # get pair counts
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF funcs
    (par, log, logf, logff, runspyder, t0) = initialize(
        "rppiC", par, nthreads=nthreads, write=write, plot=plot
    )

    # Find number of particles in input table and set nr of threads  ----------
    npt, npt1 = len(tab), len(tab1)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1])

    # Create bins in p(projected) and v(line-of-sight) space ------------------
    sepp, seppout = makebins(par.nsepp, par.seppmin, par.dsepp, par.logsepp)
    sepv, sepvout = makebins(par.nsepv, 0.0, par.dsepv, False)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei
    cra1, cdec1, cred1, cwei1 = par.cra1, par.cdec1, par.cred1, par.cwei1

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc1 = comdis(tab1[cred1].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab1 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data
        dc1 = tab1[par.cdcom1].data

    # Define the boundaries of the samples  -----------------------------------
    par.sbound = bound3d([tab[cdec].data, tab1[cdec1].data], [dc, dc1])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0".ljust(lj) + " : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # If requested, try to find the best skip grid size  ----------------------
    if par.autogrid:
        log.info("SK Autogrid".ljust(lj) + " : ON")
        par.mxh1, par.mxh2, par.mxh3, tdens = bestSKgrid3d(
            par, [npt, npt1], [tab[cra].data, tab1[cra1].data], dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("SK Autogrid".ljust(lj) + " : OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par.mxh1, par.mxh2, par.mxh3]))

    # Sort data table/s according to some order  ------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par)
    tab = tab[sidx]
    dc = dc[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1, cred1], par)
    tab1 = tab1[sidx1]
    dc1 = dc1[sidx1]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -------------------------------------------
    tstart = time.time()
    sk1, ll1 = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt1, tab1[cra1], tab1[cdec1], dc1, par.sbound, sepv, par.nsepv
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit_dr = wunit and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    par.cntid = "DR"
    log.info("====  Counting " + par.cntid + " pairs in " + str(par.mxh1) + " DEC bands")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt = pairs_cross(par, wunit_dr, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1, dc=dc, dc1=dc1)
    tend = time.time()
    logtimming(log, par.cntid, tend - tstart)
    tacc = tend - tstart
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts and integrate along pi direction  ---------------------
    dr, bdr = tidy_counts(tt, par)
    intpi = dr.sum(axis=1)
    intpib = bdr.sum(axis=1) if par.doboot else None

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("Counts plot")
            plotcf(seppout[1], intpi, 0, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutputC(
        par, npts=[npt, npt1], binslmr=seppout, dr=dr, bootc=bdr, intpi=intpi, intpib=intpib
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_rppicounts(*seppout, dr, par, cntid="rppi")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def pcf(tab, tab1, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given two astropy tables corresponding to **data** and **random** samples,
    calculate the **two-point projected auto-correlation function (pcf)**

    All input parameters that control binning, cosmology, estimator, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles (D)

    tab1 : astropy table
        Table with random particles (R)

    par : Munch dictionary
        Input parameters. See :ref:`indic-pcf` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for DD/RR/DR counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate the CF plot on screen (and saved to disk when ``write=True``)
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dd, counts.rr, counts.wrp, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicpcf` for a
        detailed description.


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                       # read data
        rans  = Table.read('redgal_rand.fits')                  # read randoms
        par = gun.packpars(kind='pcf',outfn='redgal')           # generate default parameters
        cnt = gun.pcf(gals, rans, par, write=True, plot=True)   # get pcf and plot
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize("pcf", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input tables and set nr of threads  ---------
    npt, npt1 = len(tab), len(tab1)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1])

    # Create bins in p(projected) and v(line-of-sight) space ------------------
    sepp, seppout = makebins(par.nsepp, par.seppmin, par.dsepp, par.logsepp)
    sepv, sepvout = makebins(par.nsepv, 0.0, par.dsepv, False)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei
    cra1, cdec1, cred1, cwei1 = par.cra1, par.cdec1, par.cred1, par.cwei1

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc1 = comdis(tab1[cred1].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab1 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data
        dc1 = tab1[par.cdcom1].data

    # Write out the boundary of the survey ------------------------------------
    par.sbound = bound3d([tab[cdec].data, tab1[cdec1].data], [dc, dc1])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # Adequate pars to DD,RR and DR counts ------------------------------------
    par_dd = deepcopy(par)
    par_dd.kind = "rppiA"
    par_dd.cntid = "DD"
    par_rr = deepcopy(par)
    par_rr.kind = "rppiA"
    par_rr.cntid = "RR"
    par_rr.wfib = False  # don't do fiber corrections for RR pairs
    par_rr.doboot = False  # don't do bootstraping for RR pairs
    par_dr = deepcopy(par)
    par_dr.kind = "rppiC"
    par_dr.cntid = "DR"
    par_dr.wfib = False  # don't do fiber corrections for RR pairs
    par_dr.doboot = False  # don't do bootstraping for RR pairs

    # If requested, try to find the best SK grid size for DD/RR counts  ------
    if par.autogrid:
        log.info("SK Autogrid".ljust(lj) + " : ON")
        par_dd.mxh1, par_dd.mxh2, par_dd.mxh3, tdens_dd = bestSKgrid3d(
            par_dd, npt, tab[cra].data, dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens_dd:0.3f}")
        par_rr.mxh1, par_rr.mxh2, par_rr.mxh3, tdens_rr = bestSKgrid3d(
            par_rr, npt1, tab1[cra1].data, dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens_rr:0.3f}")
        # For crosscounts choose the grid of randoms. Change if passing Random-Data order instead
        par_dr.mxh1, par_dr.mxh2, par_dr.mxh3 = par_rr.mxh1, par_rr.mxh2, par_rr.mxh3
        # par_dr.mxh1, par_dr.mxh2, par_dr.mxh3 = par_dd.mxh1, par_dd.mxh2, par_dd.mxh3
    else:
        log.info("SK Autogrid".ljust(lj) + " : OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par_dd.mxh1, par_dd.mxh2, par_dd.mxh3]))
    log.info("SK grid1 size [dec,ra,dcom]".ljust(lj) + " : " + str([par_rr.mxh1, par_rr.mxh2, par_rr.mxh3]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par_dd)
    tab, dc = tab[sidx], dc[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1, cred1], par_rr)
    tab1, dc1 = tab1[sidx1], dc1[sidx1]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -----------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll3d(
        par_dd.mxh1, par_dd.mxh2, par_dd.mxh3, npt, tab[cra], tab[cdec], dc, par_dd.sbound, sepv, par_dd.nsepv
    )
    sk1, ll1 = cff.mod.skll3d(
        par_rr.mxh1,
        par_rr.mxh2,
        par_rr.mxh3,
        npt1,
        tab1[cra1],
        tab1[cdec1],
        dc1,
        par_rr.sbound,
        sepv,
        par_rr.nsepv,
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------------------------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit_dr = wunit and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    log.info("====  Counting " + par_dd.cntid + " pairs in " + str(par_dd.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_auto(par_dd, wunit, logff, tab, x, y, z, sk, ll, dc=dc)
    tend = time.time()
    logtimming(log, par_dd.cntid, tend - tstart)
    tacc = tend - tstart

    if par.estimator not in ("DP"):
        log.info("====  Counting " + par_rr.cntid + " pairs in " + str(par_rr.mxh1) + " DEC strips  =====")
        if runspyder:
            log.info("      [for progress updates check " + logff + "]")
        tstart = time.time()
        tt2 = pairs_auto(par_rr, wunit1, logff, tab1, x1, y1, z1, sk1, ll1, dc=dc1)
        tend = time.time()
        logtimming(log, par_rr.cntid, tend - tstart)
        tacc = tacc + (tend - tstart)

    if par.estimator in ("HAM", "DP", "LS"):
        log.info("====  Counting " + par_dr.cntid + " pairs in " + str(par_dr.mxh1) + " DEC strips  =====")
        if runspyder:
            log.info("      [for progress updates check " + logff + "]")
        tstart = time.time()
        tt3 = pairs_cross(par_dr, wunit_dr, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1, dc=dc, dc1=dc1)
        tend = time.time()
        logtimming(log, par_dr.cntid, tend - tstart)
        tacc = tacc + (tend - tstart)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    dd, bdd = tidy_counts(tt1, par_dd)
    rr, dum = (
        tidy_counts(tt2, par_rr) if par.estimator not in ("DP") else (np.zeros([par.nsepp, par.nsepv]), 0)
    )
    dr, dum = (
        tidy_counts(tt3, par_dr)
        if par.estimator in ("HAM", "DP", "LS")
        else (np.zeros([par.nsepp, par.nsepv]), 0)
    )

    # Compute projected correlation function estimate  ------------------------
    (wrp, wrperr) = tpcf_wrp(npt, npt1, dd, bdd, rr, dr, par.dsepv, estimator=par.estimator)
    # rpl, rpm, rpr = seppout

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("PCF plot 1")
            plotcf(seppout[1], wrp, wrperr, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutput(
        par, npts=[npt, npt1], binslmr=seppout, dd=dd, rr=rr, dr=dr, bootc=bdd, cf=wrp, cferr=wrperr
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_cf(*seppout, wrp, wrperr, par)
        writeasc_rppicounts(*seppout, dd, par, cntid="dd")
        if par.estimator not in ("DP"):
            writeasc_rppicounts(*seppout, rr, par, cntid="rr")
        if par.estimator in ("HAM", "DP", "LS"):
            writeasc_rppicounts(*seppout, dr, par, cntid="dr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def pccf(tab, tab1, tab2, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given three astropy tables corresponding to **data**, **random**, and
    **cross** samples, this routine calculates the **two-point projected
    cross-correlation function (pccf)**

    All input parameters that control binning, cosmology, estimator, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles (D)

    tab1 : astropy table
        Table with random particles (R)

    tab2 : astropy table
        Table with cross sample particles (C)

    par : Munch dictionary
        Input parameters. See :ref:`indic-pccf` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for CD/CR counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate the CF plot on screen (and saved to disk when ``write=True``)
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.cd, counts.cr, counts.wrp, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicpccf` for a
        detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                       # read data files
        rans  = Table.read('redgal_rand.fits')
        qsos  = Table.read('qso.fits')
        par = gun.packpars(kind='pccf',outfn='redgal_qso')      # generate default parameters
        cnt = gun.pccf(gals, rans, qsos, par, write=True )      # get pcf
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize(
        "pccf", par, nthreads=nthreads, write=write, plot=plot
    )

    # Find number of particles in input tables and set nr of threads  ---------
    npt, npt1, npt2 = len(tab), len(tab1), len(tab2)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1, npt2])

    # Create bins in p(projected) and v(line-of-sight) space ------------------
    sepp, seppout = makebins(par.nsepp, par.seppmin, par.dsepp, par.logsepp)
    sepv, sepvout = makebins(par.nsepv, 0.0, par.dsepv, False)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei
    cra1, cdec1, cred1, cwei1 = par.cra1, par.cdec1, par.cred1, par.cwei1
    cra2, cdec2, cred2, cwei2 = par.cra2, par.cdec2, par.cred2, par.cwei2

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc1 = comdis(tab1[cred1].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab1 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc2 = comdis(tab2[cred2].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab2 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data
        dc1 = tab1[par.cdcom1].data
        dc2 = tab2[par.cdcom2].data

    # Write out the boundary of the survey ------------------------------------
    par.sbound = bound3d([tab[cdec].data, tab1[cdec1].data, tab2[cdec2].data], [dc, dc1, dc2])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # Adequate pars to CD and CR counts ---------------------------------------
    par_cd = deepcopy(par)
    par_cd.kind = "rppiC"
    par_cd.cntid = "CD"
    par_cr = deepcopy(par)
    par_cr.kind = "rppiC"
    par_cr.cntid = "CR"
    par_cr.wfib = False  # don't do fiber corrections in crounts counts ?
    par_cr.doboot = False  # don't do bootstraping in cross counts ?

    # If requested, try to find the best skip grid size  ----------------------
    if par.autogrid:
        log.info("SK Autogrid ON")
        # We choose to use a single SK grid for all tables, instead of two different
        # for cd and cr counts. By feeding bestSKgrid3d() with the combination of
        # the 2 largest samples among C,D,R, we get a very good (h1,h2,h3) set.
        # Most likely, its dominated by the random sample R when npt1 is large
        ltn = [npt, npt1, npt2]
        mxposA = np.argmax(ltn)
        ltn[mxposA] = -1
        mxposB = np.argmax(ltn)
        ltn = [npt, npt1, npt2]
        ltras = [tab[cra].data, tab1[cra1].data, tab2[cra2].data]
        argn = [ltn[i] for i in [mxposA, mxposB]]
        argras = [ltras[i] for i in [mxposA, mxposB]]
        par.mxh1, par.mxh2, par.mxh3, tdens = bestSKgrid3d(par_cr, argn, argras, dens=par.dens)
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
        par_cd.mxh1, par_cd.mxh2, par_cd.mxh3 = par.mxh1, par.mxh2, par.mxh3
        par_cr.mxh1, par_cr.mxh2, par_cr.mxh3 = par.mxh1, par.mxh2, par.mxh3
    else:
        log.info("Autogrid OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par.mxh1, par.mxh2, par.mxh3]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par)  # par or par_cd, par_cr ??? XXXXX
    tab, dc = tab[sidx], dc[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1, cred1], par)
    tab1, dc1 = tab1[sidx1], dc1[sidx1]
    sidx2 = pixsort(tab2, [cra2, cdec2, cred2], par)
    tab2, dc2 = tab2[sidx2], dc2[sidx2]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -----------------------------------------------
    tstart = time.time()  # again par or par_cd ?????? XXXXXX
    sk, ll = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt, tab[cra], tab[cdec], dc, par.sbound, sepv, par.nsepv
    )
    sk1, ll1 = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt1, tab1[cra1], tab1[cdec1], dc1, par.sbound, sepv, par.nsepv
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)
    x2, y2, z2 = radec2xyz(tab2[cra2].data * np.pi / 180.0, tab2[cdec2].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit2 = (tab2[cwei2].data == 1).all()
    wunit_cd = wunit2 and wunit
    wunit_cr = wunit2 and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    log.info("====  Counting " + par_cd.cntid + " pairs in " + str(par_cd.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_cross(par_cd, wunit_cd, logff, tab2, x2, y2, z2, tab, x, y, z, sk, ll, dc=dc2, dc1=dc)
    tend = time.time()
    logtimming(log, par_cd.cntid, tend - tstart)
    tacc = tend - tstart

    log.info("====  Counting " + par_cr.cntid + " pairs in " + str(par_cr.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt2 = pairs_cross(par_cr, wunit_cr, logff, tab2, x2, y2, z2, tab1, x1, y1, z1, sk1, ll1, dc=dc2, dc1=dc1)
    tend = time.time()
    logtimming(log, par_cd.cntid, tend - tstart)
    tacc = tacc + (tend - tstart)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    cd, bcd = tidy_counts(tt1, par_cd)
    cr, dum = tidy_counts(tt2, par_cr)

    # Compute projected correlation function estimate  ------------------------
    (wrp, wrperr) = tpccf_wrp(npt, npt1, cd, bcd, cr, par.dsepv, estimator=par.estimator)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("PCCF plot 1")
            plotcf(seppout[1], wrp, wrperr, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput  ------------------------------------------------------------
    counts = buildoutput(
        par, npts=[npt, npt1, npt2], binslmr=seppout, cd=cd, cr=cr, bootc=bcd, cf=wrp, cferr=wrperr
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_cf(*seppout, wrp, wrperr, par)
        writeasc_rppicounts(*seppout, cd, par, cntid="cd")
        writeasc_rppicounts(*seppout, cr, par, cntid="cr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def s_A(tab, par, nthreads=-1, write=True, para=False, plot=False, **kwargs):
    """
    Given an astropy table, count pairs in redshift space

    All input parameters that control binning, cosmology, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles

    par : Munch dictionary
        Input parameters. See :ref:`indic-sA` for a detailed description

    write : bool. Default=True
       * True : write several output files for counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    para : bool. Default=False
       * True : run in parallel over all available ipyparallel engines
       * False : run serially. No need for any ipyparallel cluster running

    plot : bool. Default=False
       * True : generate a plot of counts vs redshift separation
       * False : do not generate any plots


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dd, counts.bdd, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicsA`
        for a detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                   # read data
        par = gun.packpars(kind='sA',outfn='redgalpairs')   # generate default parameters
        cnt = gun.s_A(gals, par, write=True, plot=True  )   # get counts and plot
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize("sA", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input table and set nr of threads  ----------
    npt = len(tab)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt])

    # Create bins in redshift space -------------------------------------------
    seps, sepsout = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data

    # Define the boundaries of the sample  ------------------------------------
    par.sbound = bound3d(tab[cdec].data, dc)
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0".ljust(lj) + " : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # If requested, try to find the best SK grid size  ------------------------
    if par.autogrid:
        log.info("SK Autogrid".ljust(lj) + " : ON")
        par.mxh1, par.mxh2, par.mxh3, tdens = bestSKgrid3d(par, npt, tab[cra].data, dens=par.dens)
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("SK Autogrid".ljust(lj) + " : OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par.mxh1, par.mxh2, par.mxh3]))

    # Sort data table/s according to some order  ------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par)
    tab, dc = tab[sidx], dc[sidx]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  ------------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt, tab[cra], tab[cdec], dc, par.sbound, seps, par.nseps
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------------------------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1.0).all()

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    par.cntid = "DD"
    log.info("====  Counting " + par.cntid + " pairs in " + str(par.mxh1) + " DEC bands  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_auto(par, wunit, logff, tab, x, y, z, sk, ll, dc=dc)
    tend = time.time()
    logtimming(log, par.cntid, tend - tstart)
    tacc = tend - tstart
    #        dv.push(dict(npt=npt, ra=tab[cra].data, dec=tab[cdec].data, dc=dc, wei=tab[cwei].data,
    #                     x=x, y=y, z=z, seps=seps, sk=sk, ll=ll, allw=allw, par=par), block=True)
    #        tt1 = fcall_sA_serial(npt, tab[cra].data, tab[cdec].data, dc, tab[cwei].data,
    #                              x, y, z, seps, sk, ll, allw, par)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy output counts  -----------------------------------------------------
    dd, bdd = tidy_counts(tt1, par)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("Counts plot")
            plotcf(sepsout[1], dd, 0, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutputC(par, npts=[npt], binslmr=sepsout, dd=dd, bootc=bdd)

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_counts(*sepsout, dd, par, cntid="dd")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def s_C(tab, tab1, par, nthreads=-1, write=True, para=False, plot=False, **kwargs):
    """
    Given two astropy tables, cross-count pairs in redshift space

    All input parameters that control binning, cosmology, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with particles (sample D)

    tab1 : astropy table
        Table with particles (sample R)

    par : Munch dictionary
        Input parameters. See :ref:`indic-sC` for a detailed description

    write : bool. Default=True
       * True : write several output files for counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    para : bool. Default=False
       * True : run in parallel over all available ipyparallel engines
       * False : run serially. No need for any ipyparallel cluster running

    plot : bool. Default=False
       * True : generate a plot of counts vs redshift separation
       * False : do not generate any plots


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dr, counts.bdr, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicsC`
        for a detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        qsos = Table.read('qso.fits')                        # read data
        gals = Table.read('redgal.fits')                     # read data
        par = gun.packpars(kind='sC',outfn='qso_rg_pairs')   # generate default parameters
        cnt = gun.rppiC(qsos, gals, par, write=True)         # get pair counts
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    t0 = time.time()
    (par, log, logf, logff, runspyder, t0) = initialize("sC", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input table and set nr of threads  ----------
    npt, npt1 = len(tab), len(tab1)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1])

    # Create bins in redshift space -------------------------------------------
    seps, sepsout = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei
    cra1, cdec1, cred1, cwei1 = par.cra1, par.cdec1, par.cred1, par.cwei1

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc1 = comdis(tab1[cred1].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab1 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data
        dc1 = tab1[par.cdcom1].data

    # Write out the boundary of the survey -------------
    par.sbound = bound3d([tab[cdec].data, tab1[cdec1].data], [dc, dc1])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # If requested, try to find the best SK grid size  ------------------------
    if par.autogrid:
        log.info("SK Autogrid".ljust(lj) + " : ON")
        par.mxh1, par.mxh2, par.mxh3, tdens = bestSKgrid3d(
            par, [npt, npt1], [tab[cra].data, tab1[cra1].data], dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("SK Autogrid".ljust(lj) + " : OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par.mxh1, par.mxh2, par.mxh3]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par)
    tab, dc = tab[sidx], dc[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1, cred1], par)
    tab1, dc1 = tab1[sidx1], dc1[sidx1]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -------------------------------------------
    tstart = time.time()
    sk1, ll1 = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt1, tab1[cra1], tab1[cdec1], dc1, par.sbound, seps, par.nseps
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit_dr = wunit and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    par.cntid = "DR"
    log.info("====  Counting " + par.cntid + " pairs in " + str(par.mxh1) + " DEC bands  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_cross(par, wunit_dr, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1, dc=dc, dc1=dc1)
    tend = time.time()
    logtimming(log, par.cntid, tend - tstart)
    tacc = tend - tstart

    #    dv.push(dict(npt=npt, ra=tab[cra].data, dec=tab[cdec].data, dc=dc, wei=tab[cwei].data,  x=x, y=y, z=z,
    #                 npt1=npt1, ra1=tab1[cra1].data, dec1=tab1[cdec1].data, dc1=dc1, wei1=tab1[cwei1].data,  x1=x1, y1=y1, z1=z1,
    #                 seps=seps,sk1=sk1, ll1=ll1, allw=allw_dr, par=par), block=True)
    #
    #
    #    tt1 = fcall_sC_serial(npt, tab[cra].data, tab[cdec].data, dc, tab[cwei].data, x, y, z,
    #                          npt1, tab1[cra1], tab1[cdec1].data, dc1, tab1[cwei1].data, x1, y1, z1,
    #                          seps, sk1, ll1, allw_dr, par)

    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    dr, bdr = tidy_counts(tt1, par)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("Counts plot")
            plotcf(sepsout[1], dr, 0, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutputC(par, npts=[npt, npt1], binslmr=sepsout, dr=dr, bootc=bdr)

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_counts(*sepsout, dr, par, cntid="dr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def rcf(tab, tab1, par, nthreads=-1, write=True, para=False, plot=False, **kwargs):
    """
    Given two astropy tables corresponding to **data** and **random** samples,
    this routine calculates the **two-point redshift space auto-correlation
    function (rcf)**

    All input parameters that control binning, cosmology, estimator, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles (D)

    tab1 : astropy table
        Table with random particles (R)

    par : Munch dictionary
        Input parameters. See :ref:`indic-rcf` for a detailed description

    write : bool. Default=True
       * True : write several output files for DD/RR/DR counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    para : bool. Default=False
       * True : run in parallel over all available ipyparallel engines
       * False : run serially. No need for any ipyparallel cluster running

    plot : bool. Default=False
       * True : generate the CF plot on screen (and saved to disk when ``write=True``)
       * False : do not generate any plots


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dd, counts.rr, counts.xis, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicrcf` for a
        detailed description.


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                       # read data
        rans  = Table.read('redgal_rand.fits')                  # read randoms
        par = gun.packpars(kind='rcf',outfn='redgal')           # generate default parameters
        cnt = gun.rcf(gals, rans, par, write=True, plot=True)   # get rcf and plot
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    t0 = time.time()
    (par, log, logf, logff, runspyder, t0) = initialize("rcf", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input tables and set nr of threads  ---------
    npt, npt1 = len(tab), len(tab1)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1])

    # Create bins in redshift space -------------------------------------------
    seps, sepsout = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei
    cra1, cdec1, cred1, cwei1 = par.cra1, par.cdec1, par.cred1, par.cwei1

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc1 = comdis(tab1[cred1].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab1 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data
        dc1 = tab1[par.cdcom1].data

    # Write out the boundary of the survey ------------------------------------
    par.sbound = bound3d([tab[cdec].data, tab1[cdec1].data], [dc, dc1])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # Adequate pars to DD,RR and DR counts ------------------------------------
    par_dd = deepcopy(par)
    par_dd.kind = "sA"
    par_dd.cntid = "DD"
    par_rr = deepcopy(par)
    par_rr.kind = "sA"
    par_rr.cntid = "RR"
    par_rr.wfib = False  # don't do fiber corrections in random counts
    par_rr.doboot = False  # don't do bootstraping in random counts
    par_dr = deepcopy(par)
    par_dr.kind = "sC"
    par_dr.cntid = "DR"
    par_dr.wfib = False  # don't do fiber corrections in crounts counts
    par_dr.doboot = False  # don't do bootstraping in cross counts

    # If requested, try to find the best SK grid for DD/RR counts  -----------
    if par.autogrid:
        log.info("SK Autogrid ON")
        par_dd.mxh1, par_dd.mxh2, par_dd.mxh3, tdens_dd = bestSKgrid3d(
            par_dd, npt, tab[cra].data, dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens_dd:0.3f}")
        par_rr.mxh1, par_rr.mxh2, par_rr.mxh3, tdens_rr = bestSKgrid3d(
            par_rr, npt1, tab1[cra1].data, dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens_rr:0.3f}")
        # For crosscounts choose the grid of randoms. Change if passing Random-Data order instead
        par_dr.mxh1, par_dr.mxh2, par_dr.mxh3 = par_rr.mxh1, par_rr.mxh2, par_rr.mxh3
        # par_dr.mxh1,par_dr.mxh2,par_dr.mxh3 = par_dd.mxh1,par_dd.mxh2,par_dd.mxh3
    else:
        log.info("SK Autogrid".ljust(lj) + " : OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par_dd.mxh1, par_dd.mxh2, par_dd.mxh3]))
    log.info("SK grid1 size [dec,ra,dcom]".ljust(lj) + " : " + str([par_rr.mxh1, par_rr.mxh2, par_rr.mxh3]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par_dd)
    tab, dc = tab[sidx], dc[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1, cred1], par_rr)
    tab1, dc1 = tab1[sidx1], dc1[sidx1]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -----------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll3d(
        par_dd.mxh1, par_dd.mxh2, par_dd.mxh3, npt, tab[cra], tab[cdec], dc, par_dd.sbound, seps, par_dd.nseps
    )
    sk1, ll1 = cff.mod.skll3d(
        par_rr.mxh1,
        par_rr.mxh2,
        par_rr.mxh3,
        npt1,
        tab1[cra1],
        tab1[cdec1],
        dc1,
        par_rr.sbound,
        seps,
        par_rr.nseps,
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------------------------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit_dr = wunit and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    log.info("====  Counting " + par_dd.cntid + " pairs in " + str(par_dd.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_auto(par_dd, wunit, logff, tab, x, y, z, sk, ll, dc=dc)
    tend = time.time()
    logtimming(log, par_dd.cntid, tend - tstart)
    tacc = tend - tstart

    if par.estimator not in ("DP"):
        log.info("====  Counting " + par_rr.cntid + " pairs in " + str(par_rr.mxh1) + " DEC strips  =====")
        if runspyder:
            log.info("      [for progress updates check " + logff + "]")
        tstart = time.time()
        tt2 = pairs_auto(par_rr, wunit1, logff, tab1, x1, y1, z1, sk1, ll1, dc=dc1)
        tend = time.time()
        logtimming(log, par_rr.cntid, tend - tstart)
        tacc = tacc + (tend - tstart)

    if par.estimator in ("HAM", "DP", "LS"):
        log.info("====  Counting " + par_dr.cntid + " pairs in " + str(par_dr.mxh1) + " DEC strips  =====")
        if runspyder:
            log.info("      [for progress updates check " + logff + "]")
        tstart = time.time()
        tt3 = pairs_cross(par_dr, wunit_dr, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1, dc=dc, dc1=dc1)
        tend = time.time()
        logtimming(log, par_dr.cntid, tend - tstart)
        tacc = tacc + (tend - tstart)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    #    tt1 = fcall_sA_serial(npt, tab[cra].data, tab[cdec].data, dc, tab[cwei].data,
    #                          x, y, z, seps, sk, ll, allw, par_dd)
    #    tt2 = fcall_sA_serial(npt1, tab1[cra1].data, tab1[cdec1].data, dc1, tab1[cwei1].data,
    #                          x1, y1, z1, seps, sk1, ll1, allw1, par_rr)
    #    tt3 = fcall_sC_serial(npt1, tab1[cra1].data, tab1[cdec1].data, dc1, tab1[cwei1].data, x1, y1, z1,
    #                          npt, tab[cra], tab[cdec].data, dc, tab[cwei].data, x, y, z,
    #                          seps, sk, ll, allw_dr, par_dr)

    # Tidy ouput counts  ------------------------------------------------------
    dd, bdd = tidy_counts(tt1, par_dd)
    rr, dum = tidy_counts(tt2, par_rr) if par.estimator not in ("DP") else (np.zeros(par.nseps), 0)
    dr, dum = tidy_counts(tt3, par_dr) if par.estimator in ("HAM", "DP", "LS") else (np.zeros(par.nseps), 0)

    # Compute projected correlation function estimate  ------------------------
    (xis, xiserr) = tpcf(npt, npt1, dd, bdd, rr, dr, estimator=par.estimator)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("RCF plot 1")
            plotcf(sepsout[1], xis, xiserr, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutput(
        par, npts=[npt, npt1], binslmr=sepsout, dd=dd, rr=rr, dr=dr, bootc=bdd, cf=xis, cferr=xiserr
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_cf(*sepsout, xis, xiserr, par)
        writeasc_counts(*sepsout, dd, par, cntid="dd")
        if par.estimator not in ("DP"):
            writeasc_counts(*sepsout, rr, par, cntid="rr")
        if par.estimator in ("HAM", "DP", "LS"):
            writeasc_counts(*sepsout, dr, par, cntid="dr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def rccf(tab, tab1, tab2, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given three astropy tables corresponding to **data**, **random**, and
    **cross** samples, this routine calculates the **two-point redshift space
    cross-correlation function (rccf)**

    All input parameters that control binning, cosmology, estimator, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles (D)

    tab1 : astropy table
        Table with random particles (R)

    tab2 : astropy table
        Table with cross sample particles (C)

    par : Munch dictionary
        Input parameters. See :ref:`indic-rccf` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for CD/CR counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate the CF plot on screen (and saved to disk when ``write=True``)
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.cd, counts.cr, counts.xis, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicrccf` for a
        detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                       # read data files
        rans  = Table.read('redgal_rand.fits')
        qsos  = Table.read('qso.fits')
        par = gun.packpars(kind='pccf',outfn='redgal_qso')      # generate default parameters
        cnt = gun.rccf(gals, rans, qsos, par, write=True )      # get rcf
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    t0 = time.time()
    (par, log, logf, logff, runspyder, t0) = initialize(
        "rccf", par, nthreads=nthreads, write=write, plot=plot
    )

    # Find number of particles in input tables and set nr of threads  ---------
    npt, npt1, npt2 = len(tab), len(tab1), len(tab2)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1, npt2])

    # Create bins in redshift space -------------------------------------------
    seps, sepsout = makebins(par.nseps, par.sepsmin, par.dseps, par.logseps)

    # Unpack column names in params just for shorter writing  ------
    cra, cdec, cred, cwei = par.cra, par.cdec, par.cred, par.cwei
    cra1, cdec1, cred1, cwei1 = par.cra1, par.cdec1, par.cred1, par.cwei1
    cra2, cdec2, cred2, cwei2 = par.cra2, par.cdec2, par.cred2, par.cwei2

    # Get comoving distances --------------------------------------------------
    if par.calcdist:
        tstart = time.time()
        dc = comdis(tab[cred].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc1 = comdis(tab1[cred1].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab1 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
        tstart = time.time()
        dc2 = comdis(tab2[cred2].data, par, nt)
        tend = time.time()
        log.info("Comov_dist_tab2 compute time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")
    else:
        log.info("Using input comov. distances")
        dc = tab[par.cdcom].data
        dc1 = tab1[par.cdcom1].data
        dc2 = tab2[par.cdcom2].data

    # Write out the boundary of the survey ------------------------------------
    par.sbound = bound3d([tab[cdec].data, tab1[cdec1].data, tab2[cdec2].data], [dc, dc1, dc2])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 6).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # Adequate pars to CD and CR counts  --------------------------------------
    par_cd = deepcopy(par)
    par_cd.kind = "sC"
    par_cd.cntid = "CD"
    par_cr = deepcopy(par)
    par_cr.kind = "sC"
    par_cr.cntid = "CR"
    par_cr.wfib = False  # don't do fiber corrections in crounts counts ?
    par_cr.doboot = False  # don't do bootstraping in cross counts ?

    # If requested, try to find the best skip grid size  ----------------------
    if par.autogrid:
        log.info("SK Autogrid ON")
        # We choose to use a single SK grid for all tables, instead of two different
        # for cd and cr counts. By feeding bestSKgrid3d() with the combination of
        # the 2 largest samples among C,D,R, we get a very good (h1,h2,h3) set.
        # Most likely, its dominated by the random sample R when npt1 is large
        ltn = [npt, npt1, npt2]
        mxposA = np.argmax(ltn)
        ltn[mxposA] = -1
        mxposB = np.argmax(ltn)
        ltn = [npt, npt1, npt2]
        ltras = [tab[cra].data, tab1[cra1].data, tab2[cra2].data]
        argn = [ltn[i] for i in [mxposA, mxposB]]
        argras = [ltras[i] for i in [mxposA, mxposB]]
        par.mxh1, par.mxh2, par.mxh3, tdens = bestSKgrid3d(par_cr, argn, argras, dens=par.dens)
        par_cd.mxh1, par_cd.mxh2, par_cd.mxh3 = par.mxh1, par.mxh2, par.mxh3
        par_cr.mxh1, par_cr.mxh2, par_cr.mxh3 = par.mxh1, par.mxh2, par.mxh3
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("Autogrid OFF")
    log.info("SK grid size [dec,ra,dcom]".ljust(lj) + " : " + str([par.mxh1, par.mxh2, par.mxh3]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec, cred], par)  # par or par_cd, par_cr ??? XXXXX
    tab, dc = tab[sidx], dc[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1, cred1], par)
    tab1, dc1 = tab1[sidx1], dc1[sidx1]
    sidx2 = pixsort(tab2, [cra2, cdec2, cred2], par)
    tab2, dc2 = tab2[sidx2], dc2[sidx2]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -----------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt, tab[cra], tab[cdec], dc, par.sbound, seps, par.nseps
    )
    sk1, ll1 = cff.mod.skll3d(
        par.mxh1, par.mxh2, par.mxh3, npt1, tab1[cra1], tab1[cdec1], dc1, par.sbound, seps, par.nseps
    )
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)
    x2, y2, z2 = radec2xyz(tab2[cra2].data * np.pi / 180.0, tab2[cdec2].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit2 = (tab2[cwei2].data == 1).all()
    wunit_cd = wunit2 and wunit
    wunit_cr = wunit2 and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    log.info("====  Counting " + par_cd.cntid + " pairs in " + str(par_cd.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_cross(par_cd, wunit_cd, logff, tab2, x2, y2, z2, tab, x, y, z, sk, ll, dc=dc2, dc1=dc)
    tend = time.time()
    logtimming(log, par_cd.cntid, tend - tstart)
    tacc = tend - tstart

    log.info("====  Counting " + par_cr.cntid + " pairs in " + str(par_cr.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt2 = pairs_cross(par_cr, wunit_cr, logff, tab2, x2, y2, z2, tab1, x1, y1, z1, sk1, ll1, dc=dc2, dc1=dc1)
    tend = time.time()
    logtimming(log, par_cd.cntid, tend - tstart)
    tacc = tacc + (tend - tstart)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================
    #    tt1 = fcall_sC_serial(npt2, tab2[cra2].data, tab2[cdec2].data, dc2, tab2[cwei2].data, x2, y2, z2,
    #                          npt, tab[cra], tab[cdec].data, dc, tab[cwei].data, x, y, z,
    #                          seps, sk, ll, allw_cd, par_cd)
    #
    #    tt2 = fcall_sC_serial(npt2, tab2[cra2].data, tab2[cdec2].data, dc2, tab2[cwei2].data, x2, y2, z2,
    #                          npt1, tab1[cra1], tab1[cdec1].data, dc1, tab1[cwei1].data, x1, y1, z1,
    #                          seps, sk1, ll1, allw_cr, par_cr)

    # Tidy ouput counts  ------------------------------------------------------
    cd, bcd = tidy_counts(tt1, par_cd)
    cr, dum = tidy_counts(tt2, par_cr)

    # Compute projected correlation function estimate  ------------------------
    (xis, xiserr) = tpccf(npt, npt1, cd, bcd, cr, estimator=par.estimator)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("RCCF plot 1")
            plotcf(sepsout[1], xis, xiserr, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput  ------------------------------------------------------------
    counts = buildoutput(
        par, npts=[npt, npt1, npt2], binslmr=sepsout, cd=cd, cr=cr, bootc=bcd, cf=xis, cferr=xiserr
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write ascii counts, correlations and parameters  ------------------------
    if write:
        writeasc_cf(*sepsout, xis, xiserr, par)
        writeasc_counts(*sepsout, cd, par, cntid="cd")
        writeasc_counts(*sepsout, cr, par, cntid="cr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def th_A(tab, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given an astropy data table, count pairs in angular space

    All input parameters that control binning, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles

    par : Munch dictionary
        Input parameters. See :ref:`indic-thA` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate a plot of counts vs angular separation
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dd, counts.bdd, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicthA`
        for a detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                   # read data
        par = gun.packpars(kind='thA',outfn='redgalpairs')  # generate default parameters
        cnt = gun.th_A(gals, par, write=True, plot=True  )  # get counts and plot
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize("thA", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input table and set nr of threads  ----------
    npt = len(tab)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt])

    # Create bins in angular space --------------------------------------------
    sept, septout = makebins(par.nsept, par.septmin, par.dsept, par.logsept)

    # Unpack column names in params just for shorter writing  ------
    cra, cdec, cwei = par.cra, par.cdec, par.cwei

    # Define the boundaries of the sample  ------------------------------------
    par.sbound = bound2d(tab[cdec].data)
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 4).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # If requested, try to find the best SK grid size  ------------------------
    if par.autogrid:
        log.info("SK Autogrid ON")
        par.mxh1, par.mxh2, tdens = bestSKgrid2d(par, npt, tab[cra].data, dens=par.dens)
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("Autogrid OFF")
    log.info("SK grid size [dec,ra]".ljust(lj) + " : " + str([par.mxh1, par.mxh2]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec], par)
    tab = tab[sidx]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll2d(par.mxh1, par.mxh2, npt, tab[cra], tab[cdec], par.sbound)
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    par.cntid = "DD"
    log.info("====  Counting " + par.cntid + " pairs in " + str(par.mxh1) + " DEC bands  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt = pairs_auto(par, wunit, logff, tab, x, y, z, sk, ll)
    tend = time.time()
    logtimming(log, par.cntid, tend - tstart)
    tacc = tend - tstart
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    dd, bdd = tidy_counts(tt, par)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("Counts plot")
            plotcf(septout[1], dd, 0, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutputC(par, npts=[npt], binslmr=septout, dd=dd, bootc=bdd)

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_counts(*septout, dd, par, cntid="dd")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def th_C(tab, tab1, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given two astropy data tables, cross-count pairs in angular space

    All input parameters that control binning, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with particles (sample D)

    tab1 : astropy table
        Table with particles (sample R)

    par : Munch dictionary
        Input parameters. See :ref:`indic-thC` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate a plot of counts vs angular separation
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dr, counts.bdr, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicthC`
        for a detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        qsos = Table.read('qso.fits')                        # read data
        gals = Table.read('redgal.fits')                     # read data
        par = gun.packpars(kind='thC',outfn='qso_rg_pairs')  # generate default parameters
        cnt = gun.th_C(qsos, gals, par, write=True)          # get pair counts
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize("thC", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input table and set nr of threads  ----------
    npt, npt1 = len(tab), len(tab1)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1])

    # Create bins in angular space --------------------------------------------
    sept, septout = makebins(par.nsept, par.septmin, par.dsept, par.logsept)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cwei = par.cra, par.cdec, par.cwei
    cra1, cdec1, cwei1 = par.cra1, par.cdec1, par.cwei1

    # Write out the boundary of the survey -------------
    par.sbound = bound2d([tab[cdec].data, tab1[cdec1].data])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 4).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # If requested, try to find the best SK grid size  ------------------------
    if par.autogrid:
        log.info("SK Autogrid ON")
        par.mxh1, par.mxh2, tdens = bestSKgrid2d(
            par, [npt, npt1], [tab[cra].data, tab1[cra1].data], dens=par.dens
        )
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("Autogrid OFF")
    log.info("SK grid size [dec,ra]".ljust(lj) + " : " + str([par.mxh1, par.mxh2]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec], par)
    tab = tab[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1], par)
    tab1 = tab1[sidx1]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -------------------------------------------
    tstart = time.time()
    sk1, ll1 = cff.mod.skll2d(par.mxh1, par.mxh2, npt1, tab1[cra1], tab1[cdec1], par.sbound)
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit_dr = wunit and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    par.cntid = "DR"
    log.info("====  Counting " + par.cntid + " pairs in " + str(par.mxh1) + " DEC bands")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt = pairs_cross(par, wunit_dr, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1)
    tend = time.time()
    logtimming(log, par.cntid, tend - tstart)
    tacc = tend - tstart
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    dr, bdr = tidy_counts(tt, par)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("Counts plot")
            plotcf(septout[1], dr, 0, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput COUNTS object  ----------------------------------------------
    counts = buildoutputC(par, npts=[npt, npt1], binslmr=septout, dr=dr, bootc=bdr)

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_counts(*septout, dr, par, cntid="dr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def acf(tab, tab1, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given two astropy tables corresponding to **data** and **random** samples,
    this routine calculates the **two-point angular space auto-correlation
    function (acf)**

    All input parameters that control binning, cosmology, estimator, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles (D)

    tab1 : astropy table
        Table with random particles (R)

    par : Munch dictionary
        Input parameters. See :ref:`indic-acf` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for DD/RR/DR counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate the CF plot on screen (and saved to disk when ``write=True``)
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.dd, counts.rr, counts.xis, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicacf` for a
        detailed description.


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                       # read data
        rans  = Table.read('redgal_rand.fits')                  # read randoms
        par = gun.packpars(kind='rcf',outfn='redgal')           # generate default parameters
        cnt = gun.rcf(gals, rans, par, write=True, plot=True)   # get rcf and plot
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize("acf", par, nthreads=nthreads, write=write, plot=plot)

    # Find number of particles in input tables --------------------------------
    npt, npt1 = len(tab), len(tab1)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1])

    # Create bins in angular space  -------------------------------------------
    sept, septout = makebins(par.nsept, par.septmin, par.dsept, par.logsept)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cwei = par.cra, par.cdec, par.cwei
    cra1, cdec1, cwei1 = par.cra1, par.cdec1, par.cwei1

    # Write out the boundary of the survey  -----------------------------------
    par.sbound = bound2d([tab[cdec].data, tab1[cdec1].data])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 4).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # Adequate pars to DD,RR and DR counts  -----------------------------------
    par_dd = deepcopy(par)
    par_dd.kind = "thA"
    par_dd.cntid = "DD"
    par_rr = deepcopy(par)
    par_rr.kind = "thA"
    par_rr.cntid = "RR"
    par_rr.wfib = False  # don't do fiber corrections in random counts
    par_rr.doboot = False  # don't do bootstraping in random counts
    par_dr = deepcopy(par)
    par_dr.kind = "thC"
    par_dr.cntid = "DR"
    par_dr.wfib = False  # don't do fiber corrections in crounts counts
    par_dr.doboot = False  # don't do bootstraping in cross counts

    # If requested, try to find the best SK grid for DD/RR counts  -----------
    if par.autogrid:
        log.info("SK Autogrid ON")
        par_dd.mxh1, par_dd.mxh2, tdens_dd = bestSKgrid2d(par_dd, npt, tab[cra].data, dens=par.dens)
        log.info("SK cell target density".ljust(lj) + f" : {tdens_dd:0.3f}")
        par_rr.mxh1, par_rr.mxh2, tdens_rr = bestSKgrid2d(par_rr, npt1, tab1[cra1].data, dens=par.dens)
        log.info("SK cell target density".ljust(lj) + f" : {tdens_rr:0.3f}")
        # For crosscounts choose the grid of randoms. Change if passing Random-Data order instead
        par_dr.mxh1, par_dr.mxh2, par_dr.mxh3 = par_rr.mxh1, par_rr.mxh2, par_rr.mxh3
        # par_dr.mxh1, par_dr.mxh2, par_dr.mxh3 = par_dd.mxh1, par_dd.mxh2, par_dd.mxh3
    else:
        log.info("Autogrid OFF")
    log.info("SK grid size [dec,ra]".ljust(lj) + " : " + str([par_dd.mxh1, par_dd.mxh2]))
    log.info("SK grid1 size [dec,ra]".ljust(lj) + " : " + str([par_rr.mxh1, par_rr.mxh2]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec], par_dd)
    tab = tab[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1], par_rr)
    tab1 = tab1[sidx1]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -----------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll2d(par_dd.mxh1, par_dd.mxh2, npt, tab[cra], tab[cdec], par_dd.sbound)
    sk1, ll1 = cff.mod.skll2d(par_rr.mxh1, par_rr.mxh2, npt1, tab1[cra1], tab1[cdec1], par_rr.sbound)
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit_dr = wunit and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    log.info("====  Counting " + par_dd.cntid + " pairs in " + str(par_dd.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_auto(par_dd, wunit, logff, tab, x, y, z, sk, ll)
    tend = time.time()
    logtimming(log, par_dd.cntid, tend - tstart)
    tacc = tend - tstart

    if par.estimator not in ("DP"):
        log.info("====  Counting " + par_rr.cntid + " pairs in " + str(par_rr.mxh1) + " DEC strips  =====")
        if runspyder:
            log.info("      [for progress updates check " + logff + "]")
        tstart = time.time()
        tt2 = pairs_auto(par_rr, wunit1, logff, tab1, x1, y1, z1, sk1, ll1)
        tend = time.time()
        logtimming(log, par_rr.cntid, tend - tstart)
        tacc = tacc + (tend - tstart)

    if par.estimator in ("HAM", "DP", "LS"):
        log.info("====  Counting " + par_dr.cntid + " pairs in " + str(par_dr.mxh1) + " DEC strips  =====")
        if runspyder:
            log.info("      [for progress updates check " + logff + "]")
        tstart = time.time()
        tt3 = pairs_cross(par_dr, wunit_dr, logff, tab, x, y, z, tab1, x1, y1, z1, sk1, ll1)
        tend = time.time()
        logtimming(log, par_dr.cntid, tend - tstart)
        tacc = tacc + (tend - tstart)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    dd, bdd = tidy_counts(tt1, par_dd)
    rr, dum = tidy_counts(tt2, par_rr) if par.estimator not in ("DP") else (np.zeros(par.nsept), 0)
    dr, dum = tidy_counts(tt3, par_dr) if par.estimator in ("HAM", "DP", "LS") else (np.zeros(par.nsept), 0)

    # Compute angular correlation function estimate  --------------------------
    (wth, wtherr) = tpcf(npt, npt1, dd, bdd, rr, dr, estimator=par.estimator)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("ACF plot 1")
            plotcf(septout[1], wth, wtherr, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput  ------------------------------------------------------------
    counts = buildoutput(
        par, npts=[npt, npt1], binslmr=septout, dd=dd, rr=rr, dr=dr, bootc=bdd, cf=wth, cferr=wtherr
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write ascii counts, correlations and parameters  ------------------------
    if write:
        writeasc_cf(*septout, wth, wtherr, par)
        writeasc_counts(*septout, dd, par, cntid="dd")
        if par.estimator not in ("DP"):
            writeasc_counts(*septout, rr, par, cntid="rr")
        if par.estimator in ("HAM", "DP", "LS"):
            writeasc_counts(*septout, dr, par, cntid="dr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def accf(tab, tab1, tab2, par, nthreads=-1, write=True, plot=False, **kwargs):
    """
    Given three astropy tables corresponding to **data**, **random**, and
    **cross** samples, this routine calculates the **two-point angular space
    cross-correlation function (accf)**

    All input parameters that control binning, estimator, column names, etc.
    are passed in a single dictionary ``par``, which can be easily generated
    with default values using :func:`gundam.packpars` and customized later

    .. rubric :: Parameters

    tab : astropy table
        Table with data particles (D)

    tab1 : astropy table
        Table with random particles (R)

    tab2 : astropy table
        Table with cross sample particles (C)

    par : Munch dictionary
        Input parameters. See :ref:`indic-accf` for a detailed description

    nthreads : integer. Default=-1
        Number of threads to use. Default to -1, which means all available cpu cores

    write : bool. Default=True
       * True : write several output files for CD/CR counts, correlations, etc.
       * False : do not write any output files (except for logs, which are always saved)

    plot : bool. Default=False
       * True : generate the CF plot on screen (and saved to disk when ``write=True``)
       * False : do not generate any plots

    kwargs : dict
        Extra keywords passed to :func:`gundam.plotcf`


    .. rubric :: Returns

    counts : Munch dictionary
        Munch dictionary, containing all counts and correlations accesible by
        field keys, e.g. counts.cd, counts.cr, counts.wth, etc. It also stores
        the complete log and all input parameters. See :ref:`outdicaccf` for a
        detailed description


    .. rubric :: Examples

    .. code-block:: python

        import gundam as gun ; from astropy.table import Table
        gals  = Table.read('redgal.fits')                       # read data files
        rans  = Table.read('redgal_rand.fits')
        qsos  = Table.read('qso.fits')
        par = gun.packpars(kind='accf',outfn='redgal_qso')      # generate default parameters
        cnt = gun.accf(gals, rans, qsos, par, write=True )      # get accf
    """
    lj = 27  # nr of characters for left justification of some status msgs

    # Initialize logs, check if par has the right kind, etc. Common to all CF functions
    (par, log, logf, logff, runspyder, t0) = initialize(
        "accf", par, nthreads=nthreads, write=write, plot=plot
    )

    # Find number of particles in input tables and set nr of threads  ---------
    npt, npt1, npt2 = len(tab), len(tab1), len(tab2)
    nt = set_threads(nthreads)

    # Log calling information  ------------------------------------------------
    logcallinfo(log, par, npts=[npt, npt1, npt2])

    # Create bins in angular space  -------------------------------------------
    sept, septout = makebins(par.nsept, par.septmin, par.dsept, par.logsept)

    # Unpack column names in params just for shorter writing  -----------------
    cra, cdec, cwei = par.cra, par.cdec, par.cwei
    cra1, cdec1, cwei1 = par.cra1, par.cdec1, par.cwei1
    cra2, cdec2, cwei2 = par.cra2, par.cdec2, par.cwei2

    # Write out the boundary of the survey ------------------------------------
    par.sbound = bound2d([tab[cdec].data, tab1[cdec1].data, tab2[cdec2].data])
    log.info("Sample boundaries : (" + ("{:0.5f}, " * 4).format(*par.sbound)[:-2] + ")")

    # Guess if the sample cross the 360-0 deg division
    cross0 = cross0guess(tab[cra].data)
    log.info("Sample seems to cross RA=0 : " + str(cross0))
    if cross0 is True:
        log.info("Custom RA boundaries : " + str(par.custRAbound))

    # Adequate pars to CD and CR counts  --------------------------------------
    par_cd = deepcopy(par)
    par_cd.kind = "thC"
    par_cd.cntid = "CD"
    par_cr = deepcopy(par)
    par_cr.kind = "thC"
    par_cr.cntid = "CR"
    par_cr.wfib = False  # don't do fiber corrections in crounts counts ?
    par_cr.doboot = False  # don't do bootstraping in cross counts ?

    # If requested, try to find the best skip grid size  ----------------------
    if par.autogrid:
        log.info("SK Autogrid ON")
        # We choose to use a single SK grid for all tables, instead of two different
        # for cd and cr counts. By feeding bestSKgrid2d() with the combination of
        # the 2 largest samples among C,D,R, we get a very good (h1,h2,h3) set.
        # Most likely, its dominated by the random sample R when npt1 is large
        ltn = [npt, npt1, npt2]
        mxposA = np.argmax(ltn)
        ltn[mxposA] = -1
        mxposB = np.argmax(ltn)
        ltn = [npt, npt1, npt2]
        ltras = [tab[cra].data, tab1[cra1].data, tab2[cra2].data]
        argn = [ltn[i] for i in [mxposA, mxposB]]
        argras = [ltras[i] for i in [mxposA, mxposB]]
        par.mxh1, par.mxh2, tdens = bestSKgrid2d(par_cr, argn, argras, dens=par.dens)
        par_cd.mxh1, par_cd.mxh2 = par.mxh1, par.mxh2
        par_cr.mxh1, par_cr.mxh2 = par.mxh1, par.mxh2
        log.info("SK cell target density".ljust(lj) + f" : {tdens:0.3f}")
    else:
        log.info("Autogrid OFF")
    log.info("SK grid size [dec,ra]".ljust(lj) + " : " + str([par.mxh1, par.mxh2]))

    # Sort table/s according to some order  -----------------------------------
    tstart = time.time()
    sidx = pixsort(tab, [cra, cdec], par_cd)  # par or par_cd, par_cr ??? XXXXX
    tab = tab[sidx]
    sidx1 = pixsort(tab1, [cra1, cdec1], par)
    tab1 = tab1[sidx1]
    sidx2 = pixsort(tab2, [cra2, cdec2], par)
    tab2 = tab2[sidx2]
    tend = time.time()
    log.info("Pixsort time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Create SK and LL tables  -----------------------------------------------
    tstart = time.time()
    sk, ll = cff.mod.skll2d(par.mxh1, par.mxh2, npt, tab[cra], tab[cdec], par.sbound)
    sk1, ll1 = cff.mod.skll2d(par.mxh1, par.mxh2, npt1, tab1[cra1], tab1[cdec1], par.sbound)
    tend = time.time()
    log.info("SK-LL tables build time (s)".ljust(lj) + f" : {tend-tstart:0.3f}")

    # Convert ra,dec,z to spherical coords ------------
    x, y, z = radec2xyz(tab[cra].data * np.pi / 180.0, tab[cdec].data * np.pi / 180.0)
    x1, y1, z1 = radec2xyz(tab1[cra1].data * np.pi / 180.0, tab1[cdec1].data * np.pi / 180.0)
    x2, y2, z2 = radec2xyz(tab2[cra2].data * np.pi / 180.0, tab2[cdec2].data * np.pi / 180.0)

    # Find out if all weights are 1.0 to later call slighly faster functions
    wunit = (tab[cwei].data == 1).all()
    wunit1 = (tab1[cwei1].data == 1).all()
    wunit2 = (tab2[cwei2].data == 1).all()
    wunit_cd = wunit2 and wunit
    wunit_cr = wunit2 and wunit1

    # ==========================================================================
    # ==========================   COUNT PAIRS   ===============================
    log.info("====  Counting " + par_cd.cntid + " pairs in " + str(par_cd.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt1 = pairs_cross(par_cd, wunit_cd, logff, tab2, x2, y2, z2, tab, x, y, z, sk, ll)
    tend = time.time()
    logtimming(log, par_cd.cntid, tend - tstart)
    tacc = tend - tstart

    log.info("====  Counting " + par_cr.cntid + " pairs in " + str(par_cr.mxh1) + " DEC strips  =====")
    if runspyder:
        log.info("      [for progress updates check " + logff + "]")
    tstart = time.time()
    tt2 = pairs_cross(par_cr, wunit_cr, logff, tab2, x2, y2, z2, tab1, x1, y1, z1, sk1, ll1)
    tend = time.time()
    logtimming(log, par_cd.cntid, tend - tstart)
    tacc = tacc + (tend - tstart)
    # =========================  END PAIR COUNTS  ==============================
    # ==========================================================================

    # Tidy ouput counts  ------------------------------------------------------
    cd, bcd = tidy_counts(tt1, par_cd)
    cr, dum = tidy_counts(tt2, par_cr)

    # Compute projected correlation function estimate  ------------------------
    (wth, wtherr) = tpccf(npt, npt1, cd, bcd, cr, estimator=par.estimator)

    # Do plot if desired  -----------------------------------------------------
    if plot:
        try:
            # plt.figure(1)
            plt.figure("ACCF plot 1")
            plotcf(septout[1], wth, wtherr, fac=1.0, write=write, par=par, **kwargs)
        except ValueError:
            print("Warning: there is a problem with the plot !!!")

    # Build ouput  ------------------------------------------------------------
    counts = buildoutput(
        par, npts=[npt, npt1, npt2], binslmr=septout, cd=cd, cr=cr, bootc=bcd, cf=wth, cferr=wtherr
    )

    # Finalize  ---------------------------------------------------------------
    finalize(log, logf, logff, tacc, t0, counts)

    # Write output files  -----------------------------------------------------
    if write:
        writeasc_cf(*septout, wth, wtherr, par)
        writeasc_counts(*septout, cd, par, cntid="cd")
        writeasc_counts(*septout, cr, par, cntid="cr")
        savepars(par)
        savecounts(counts)

    # Close log  --------------------------------------------------------------
    closelog(log, runspyder=runspyder)

    return counts


# =============================================================================
def buildoutputC(par, npts=[], binslmr=[], dd=None, dr=None, bootc=None, intpi=None, intpib=None):
    """
    Given a set of pair count arrays, assemble the final **counts** output
    dictionary, adding also bins, input parameters and relevant sample sizes.

    This is for the main Gundam functions that calculate a single pair count, i.e.
    :func:`gundam.rppi_A`, :func:`gundam.rppi_C`, :func:`gundam.s_A`,
    :func:`gundam.s_C`, :func:`gundam.th_A`, :func:`gundam.th_C`

    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters

    npts : list of int. Default=[]
        A list with the nr. of points in each sample, i.e. [npt, npt1]. For
        example npts=[400, 4000] in cross-count cases; npts=[400] for auto-count
        cases

    binslmr : list. Default=[]
        The left, mid and right-side locations of each bin, returned of example
        by :func:`gundam.makebins`

    dd,dr : array. Default=None
        The "dd" and "dr" count arrays, if available

    bootc : array. Default=None
        The boostrap count array, if available

    intpi, intpib : array. Default=None
        "dd" counts and boostrap "dd" counts integrated along all radial bins,
        if available

    .. rubric :: Returns

    counts : Munch dictionary
        The ouput dictionary
    """

    def strz(i):
        if i == 0:
            return ""
        else:
            return str(i)

    iscross = par.kind in ["rppiC", "sC", "thC"]

    counts = Munch()
    counts.par = par  # save params dictionary
    for i in range(len(npts)):
        counts["npt" + strz(i)] = npts[i]  # save npt,npt1,npt2,...

    if par.kind in ["rppiA", "rppiC"]:  # save bins, intpi, etc.
        counts.rpl = binslmr[0]
        counts.rpm = binslmr[1]
        counts.rpr = binslmr[2]
        counts.intpi = intpi
        if par.doboot:
            counts.intpib = intpib
    elif par.kind in ["sA", "sC"]:
        counts.sl = binslmr[0]
        counts.sm = binslmr[1]
        counts.sr = binslmr[2]
    elif par.kind in ["thA", "thC"]:
        counts.thl = binslmr[0]
        counts.thm = binslmr[1]
        counts.thr = binslmr[2]
    if iscross:  # save counts
        counts.dr = dr
        if par.doboot:
            counts.bdr = bootc
    else:
        counts.dd = dd
        if par.doboot:
            counts.bdd = bootc
    return counts


# =============================================================================
def buildoutput(
    par, npts=[], binslmr=[], dd=None, rr=None, dr=None, cd=None, cr=None, bootc=None, cf=None, cferr=None
):
    """
    Given a set of pair count arrays, assemble the final **counts** output
    dictionary, adding also bins, input parameters and relevant sample sizes.

    This is for the main Gundam functions that calculate a multiple pair counts,
    i.e. :func:`gundam.pcf`, :func:`gundam.pccf`, :func:`gundam.rcf`,
    :func:`gundam.rccf`, :func:`gundam.acf`, :func:`gundam.accf`

    .. rubric :: Parameters

    par : Munch dictionary
        Input parameters

    npts : list of int. Default=[]
        A list with the nr. of points in each sample, i.e. [npt, npt1]. For
        example npts=[400, 4000] in cross-count cases; npts=[400] for auto-count
        cases

    binslmr : list. Default=[]
        The left, mid and right-side locations of each bin, returned of example
        by :func:`gundam.makebins`

    dd,rr,dr : array. Default=None
        The "dd", "rr" and "dr" count arrays, if available

    cd,cr : array. Default=None
        The "cd" and "cr" count arrays, if available

    bootc : array. Default=None
        The boostrap count array, if available

    cf,cferr : array. Default=None
        The correlation function and its error

    .. rubric :: Returns

    counts : Munch dictionary
        The ouput dictionary
    """

    def strz(i):
        if i == 0:
            return ""
        else:
            return str(i)

    iscross = par.kind in ["pccf", "rccf", "accf"]

    counts = Munch()
    counts.par = par  # save params dictionary
    for i in range(len(npts)):
        counts["npt" + strz(i)] = npts[i]  # save npt,npt1,npt2,...

    if par.kind in ["pcf", "pccf"]:  # save cf, cfeerr and the bins
        counts.wrp = cf
        counts.wrperr = cferr
        counts.rpl = binslmr[0]
        counts.rpm = binslmr[1]
        counts.rpr = binslmr[2]
    elif par.kind in ["rcf", "rccf"]:
        counts.xis = cf
        counts.xiserr = cferr
        counts.sl = binslmr[0]
        counts.sm = binslmr[1]
        counts.sr = binslmr[2]
    elif par.kind in ["acf", "accf"]:
        counts.wth = cf
        counts.wtherr = cferr
        counts.thl = binslmr[0]
        counts.thm = binslmr[1]
        counts.thr = binslmr[2]
    if iscross:
        counts.cd = cd
        counts.cr = cr
        if par.doboot:
            counts.bcd = bootc
    else:
        counts.dd = dd
        if par.estimator not in ("DP"):
            counts.rr = rr
        if par.doboot:
            counts.bdd = bootc
        if par.estimator in ("HAM", "DP", "LS"):
            counts.dr = dr

    return counts


# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
    # testsuite = 'proj_cf'                   # example pcf direct input
    testsuite = None

    if testsuite is None:
        print("Basic usage: ")
        print("       import gundam as gun")
        print("       # Read data, randoms, and create default parameters for a PCF")
        print('       gals = Table.read("galaxies.fits")')
        print('       rans = Table.read("randoms.fits")')
        print('       par  = gun.packpars(kind="pcf")')
        print("       # Calculate the correlation")
        print("       cnt = gun.pcf(gals, rans, par)")
