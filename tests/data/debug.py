import gundam as gun
from astropy.table import Table

# DEFINE PARAMETERS  ==========================================================
galf = "/export/home/phys/claudio-CMU/gundam/tests/data/DATA.fits"  # Galaxy sample
ranf = "/export/home/phys/claudio-CMU/gundam/tests/data/RAND.fits"  # Random sample

par = gun.packpars(kind="pcf")
par.autogrid = True  # Automatic SK grid size
par.mxh1 = 20  # SK size in dec
par.mxh2 = 100  # SK size in ra
par.mxh3 = 10  # SK size in z
par.nsepp = 28  # Number of bins of projected separation rp
par.seppmin = 0.02  # Minimum rp in Mpc/h
par.dsepp = 0.12  # Bin size of rp (in log space)
par.nsepv = 1  # Number of bins of LOS separation pi
par.dsepv = 40.0  # Bin size of pi (in linear space)
par.doboot = False  # Do bootstrap error estimates
par.omegam = 0.25  # Omega matter
par.omegal = 0.75  # Omega lambda
par.h0 = 100  # Hubble constant [km/s/Mpc]
par.calcdist = True  # Calculate comov. dist.
par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF


# READ DATA FILES  ============================================================
gals = Table.read(galf)
if "wei" not in gals.colnames:
    gals["wei"] = 1.0  # If not present, set weights to 1

rans = Table.read(ranf)
if "wei" not in rans.colnames:
    rans["wei"] = 1.0  # If not present, set weights to 1

# ==============================================================================
# CALCULATE THE CORRELATION
nt = 1  # Threads to use
cnt = gun.pcf(gals, rans, par, nthreads=nt, plot=False, write=False)
# ==============================================================================
