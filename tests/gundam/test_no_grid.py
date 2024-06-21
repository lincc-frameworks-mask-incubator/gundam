import gundam as gun
import numpy as np
from astropy.table import Table
from munch import Munch


def test_example_acf1():
    # DEFINE PARAMETERS  ==========================================================
    par = gun.packpars(kind="acf1")
    par.doboot = False  # Do bootstrap error estimates
    par.dsept = 0.10
    par.nsept = 33
    par.septmin = 0.01
    par.grid = 0
    par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF
    par.description = "example01"  # Description label
    par.autogrid = False

    # READ DATA FILES  ============================================================
    gals = Table.read("tests/data/DATA.fits")
    gals = gals[["ra", "dec"]]
    gals["wei"] = 1.0

    rans = Table.read("tests/data/RAND.fits")
    rans = rans[["ra", "dec"]]
    rans["wei"] = 1.0

    # ==============================================================================
    # CALCULATE THE CORRELATION
    nt = 1  # Threads to use
    cnt = gun.acf1(gals, rans, par, nthreads=nt, plot=False, write=False)
    # ==============================================================================


    r_binedges_acf = np.load("tests/data/r_binedges_acf.npy")
    c_bin_acf = np.load("tests/data/c_binedges_acf.npy")
    l_binedges_acf = np.load("tests/data/l_binedges_acf.npy")
    rr_acf = np.load("tests/data/rr_acf.npy")
    dd_acf = np.load("tests/data/dd_acf.npy")
    w_acf_nat = np.load("tests/data/w_acf_nat.npy")

    assert len(c_bin_acf) == par.nsept, f"Expected {len(c_bin_acf)} bin edges, but got {par.nsept}"
    np.testing.assert_almost_equal(
        c_bin_acf, cnt.thm, decimal=5, err_msg="Bin centers do not match expected values"
    )
    np.testing.assert_almost_equal(
        r_binedges_acf, cnt.thr, decimal=5, err_msg="Bin right borders do not match expected values"
    )
    np.testing.assert_almost_equal(
        l_binedges_acf, cnt.thl, decimal=5, err_msg="Bin left borders do not match expected values"
    )
    assert isinstance(cnt.wth, np.ndarray), f"Expected hist to be a numpy array, but got {type(np.ndarray)}"
    assert isinstance(cnt.rr, np.ndarray), f"Expected hist to be a numpy array, but got {type(np.ndarray)}"
    assert isinstance(cnt.dd, np.ndarray), f"Expected hist to be a numpy array, but got {type(np.ndarray)}"
    #assert np.issubdtype(cnt.wth.dtype, np.float), "cnt.wth dtype is not float"
    #assert np.issubdtype(cnt.rr.dtype, np.float), "cnt.rr dtype is not float"
    #assert np.issubdtype(cnt.dd.dtype, np.float), "cnt.dd dtype is not float"
    np.testing.assert_allclose(
        w_acf_nat, cnt.wth, atol=1e-1, err_msg="Correlation function is not correct"
    )
    np.testing.assert_allclose(
        rr_acf, cnt.rr, atol=2e-3, err_msg="Random-Random histogram is not correct"
    )
    np.testing.assert_allclose(
        dd_acf, cnt.dd, atol=2e-3, err_msg="Object-Object histogram is not correct"
    )

def test_example_acf_naiveway():
    # DEFINE PARAMETERS  ==========================================================
    par = Munch()
    par.dsept = 0.10
    par.nsept = 33
    par.septmin = 0.01
    par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF

    # READ DATA FILES  ============================================================
    gals = Table.read("tests/data/DATA.fits")
    gals = gals[["ra", "dec"]]

    rans = Table.read("tests/data/RAND.fits")
    rans = rans[["ra", "dec"]]
    # ==============================================================================
    # CALCULATE THE CORRELATION
    dd, rr, wth = gun.acf_naiveway(gals['ra'].data, gals['dec'].data, rans['ra'].data, rans['dec'].data, par)

    rr_acf = np.load("tests/data/rr_acf.npy")
    dd_acf = np.load("tests/data/dd_acf.npy")
    w_acf_nat = np.load("tests/data/w_acf_nat.npy")


    np.testing.assert_allclose(
        w_acf_nat, wth, atol=1e-1, err_msg="Correlation function is not correct"
    )
    np.testing.assert_allclose(
        rr_acf, rr, atol=2e-3, err_msg="Random-Random histogram is not correct"
    )
    np.testing.assert_allclose(
        dd_acf, dd, atol=2e-3, err_msg="Object-Object histogram is not correct"
    )