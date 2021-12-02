from astropy.table import Table
import gundam as gun
galf = "./data/DATA.fits"  # Galaxy sample
ranf = "./data/RAND.fits"  # Random sample
par = gun.packpars(kind="pcf", file=galf, file1=ranf, outfn='pcf_ham')
par.autogrid = False  # Automatic SK grid size
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
par.estimator = "HAM"  # Choose Landy-Szalay estimator for the PCF
#par.description = "example01"  # Description label
gals = Table.read(galf)
gals["wei"] = 1.0  # If not present, set weights to 1
rans = Table.read(ranf)
rans["wei"] = 1.0  # If not present, set weights to 1
cntnew = gun.pcf(gals, rans, par, nthreads=1, plot=False, write=True)
