from astropy.table import Table

import gundam as gun

import numpy.testing as npt

import pytest

import os
import pathlib

from pytest_easyMPI import mpi_parallel

# ============================================================================
# CONSTANTS
# ============================================================================

PATH = os.path.abspath(os.path.dirname(__file__))
RESOURCES_PATH = pathlib.Path(PATH) / "data"

@pytest.fixture
def pcf_nat():
    return gun.readcounts(RESOURCES_PATH / 'pcf_nat.cnt')

@pytest.fixture
def pcf_ham():
    return gun.readcounts(RESOURCES_PATH / 'pcf_ham.cnt')

@pytest.fixture
def pcf_ls():
    return gun.readcounts(RESOURCES_PATH / "pcf_ls.cnt")


# ============================================================================
# TESTS
# ============================================================================


# Test to prove pcf
# threads = 2
# estimator = NAT
# bootstrap = False
def test_pcf_nth2_nat(pcf_nat):
    cntoriginal = pcf_nat
    galf = str(RESOURCES_PATH) + '/DATA.fits'  # Galaxy sample
    ranf = str(RESOURCES_PATH) + '/RAND.fits'  # Random sample
    #    outfn  = './examples_data/ex_pcf_01'     # Name for output files
    par = gun.packpars(kind="pcf", file=galf, file1=ranf)
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
    par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF
    par.description = "example01"  # Description label
    gals = Table.read(galf)
    gals["wei"] = 1.0  # If not present, set weights to 1
    rans = Table.read(ranf)
    rans["wei"] = 1.0  # If not present, set weights to 1
    cntnew = gun.pcf(gals, rans, par, nthreads=2, plot=False, write=False)
    npt.assert_array_equal(cntnew.dd, cntoriginal.dd)
    npt.assert_array_equal(cntnew.rr, cntoriginal.rr)
    npt.assert_array_equal(cntnew.wrp, cntoriginal.wrp)


# Test to prove pcf
# threads = 2
# estimator = HAM
# bootstrap = False
def test_pcf_nth2_ham(pcf_ham):
    cntoriginal = pcf_ham
    galf = str(RESOURCES_PATH) + '/DATA.fits'   # Galaxy sample
    ranf = str(RESOURCES_PATH) + '/RAND.fits' # Random sample
    par = gun.packpars(kind="pcf", file=galf, file1=ranf)
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
    gals = Table.read(galf)
    gals["wei"] = 1.0  # If not present, set weights to 1
    rans = Table.read(ranf)
    rans["wei"] = 1.0  # If not present, set weights to 1
    cntnew = gun.pcf(gals, rans, par, nthreads=1, plot=False, write=False)
    npt.assert_array_equal(cntnew.dd, cntoriginal.dd)
    npt.assert_array_equal(cntnew.rr, cntoriginal.rr)
    npt.assert_array_equal(cntnew.dr, cntoriginal.dr)
    npt.assert_array_equal(cntnew.wrp, cntoriginal.wrp)

    
