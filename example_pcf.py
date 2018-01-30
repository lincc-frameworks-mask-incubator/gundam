# -*- coding: utf-8 -*-
"""
Example code to estimate a projected correlation function
"""
from astropy.table import Table
import gundam as gun

# DEFINE PARAMETERS  ==========================================================
galf   = './examples_data/DATA.fits'     # Galaxy sample
ranf   = './examples_data/RAND.fits'     # Random sample
outfn  = './examples_data/ex_pcf_01'     # Name for output files
par = gun.packpars(kind='pcf', file=galf, file1=ranf, outfn=outfn)
par.autogrid    = False       # Automatic SK grid size
par.mxh1        = 60          # SK size in dec
par.mxh2        = 240         # SK size in ra
par.mxh3        = 24          # SK size in z
par.nsepp       = 28          # Number of bins of projected separation rp
par.seppmin     = 0.02        # Minimum rp in Mpc/h
par.dsepp       = 0.12        # Bin size of rp (in log space)
par.nsepv       = 1           # Number of bins of LOS separation pi
par.dsepv       = 40.         # Bin size of pi (in linear space)
par.doboot      = True        # Do bootstrap error estimates
par.omegam      = 0.25        # Omega matter
par.omegal      = 0.75        # Omega lambda
par.h0          = 100         # Hubble constant [km/s/Mpc]
par.calcdist    = True        # Calculate comov. dist.
par.estimator   = 'LS'        # Choose Landy-Szalay estimator for the PCF
par.description = 'example01' # Description label

# READ DATA FILES  ============================================================
print('Reading file: ', galf)
gals = Table.read(galf)
if 'wei' not in gals.colnames:  gals['wei'] = 1. #If not present, set weights to 1

print('Reading file: ', ranf)
rans = Table.read(ranf)
if 'wei' not in rans.colnames:  rans['wei'] = 1. #If not present, set weights to 1

#==============================================================================
# CALCULATE THE CORRELATION
nt  = 4   # Threads to use
cnt = gun.pcf(gals, rans, par, nthreads=nt, plot=True, write=True)
#==============================================================================

