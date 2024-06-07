"""
Example code to estimate a 2D projected correlation function using a sample
of SDSS Luminous Red Galaxies
"""
from astropy.table import Table
from gundam import gundam as gun

# DEFINE PARAMETERS  ==========================================================
galf = "tests/data/DR7-lrg.fits"  # Galaxy sample
ranf = "tests/data/DR7-lrg-rand.fits"  # Random sample
outfn = "tests/data/ex_LRG"  # Name for output files


def test_example_lrg():
    par = gun.packpars(kind="pcf", file=galf, file1=ranf, outfn=outfn)
    par.autogrid = False  # Automatic SK grid size
    par.mxh1 = 60  # SK size in dec
    par.mxh2 = 240  # SK size in ra
    par.mxh3 = 24  # SK size in z
    par.nsepp = 78  # Number of bins of projected separation rp
    par.seppmin = 0.01  # Minimum rp in Mpc/h
    par.dsepp = 0.5  # Bin size of rp (in log space)
    par.logsepp = 0  # Use linear bins instead of log bins
    par.nsepv = 60  # Number of bins of LOS separation pi
    par.dsepv = 0.5  # Bin size of pi (in linear space)
    par.doboot = True  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "LS"  # Choose Landy-Szalay estimator for the PCF
    par.description = "LumRedGxs"  # Description label

    # READ DATA FILES  ============================================================
    print("Reading file: ", galf)
    gals = Table.read(galf)
    if "wei" not in gals.colnames:
        gals["wei"] = 1.0  # If not present, set weights to 1

    print("Reading file: ", ranf)
    rans = Table.read(ranf)
    if "wei" not in rans.colnames:
        rans["wei"] = 1.0  # If not present, set weights to 1

    # ==============================================================================
    # CALCULATE THE CORRELATION
    nt = 4  # Threads to use
    cnt = gun.pcf(gals, rans, par, nthreads=nt)
    # ==============================================================================

    # NOW PLOT THE 2-DIMENSIONAL PCF  =============================================
    # gun.cntplot2D(cnt, write=True)
