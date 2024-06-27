"""
Example code to estimate a projected correlation function
"""
import gundam as gun
from astropy.table import Table

# DEFINE PARAMETERS  ==========================================================
galf = "tests/data/DATA.fits"  # Galaxy sample
ranf = "tests/data/RAND.fits"  # Random sample
outfn = "tests/data/ex_pcf_01"  # Name for output files


def test_example_pcf():
    par = gun.packpars(kind="pcf", file=galf, file1=ranf, outfn=outfn)
    par.autogrid = True  # Automatic SK grid size
    par.mxh1 = 20  # SK size in dec
    par.mxh2 = 100  # SK size in ra
    par.mxh3 = 10  # SK size in z
    par.nsepp = 28  # Number of bins of projected separation rp
    par.seppmin = 0.02  # Minimum rp in Mpc/h
    par.dsepp = 0.12  # Bin size of rp (in log space)
    par.nsepv = 1  # Number of bins of LOS separation pi
    par.dsepv = 40.0  # Bin size of pi (in linear space)
    par.doboot = True  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "LS"  # Choose Landy-Szalay estimator for the PCF
    par.description = "example01"  # Description label

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
    cnt = gun.pcf(gals, rans, par, nthreads=nt, plot=True, write=True)
    print(cnt)
    # ==============================================================================
