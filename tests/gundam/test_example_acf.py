"""
Example code to estimate a projected correlation function
"""
import gundam as gun
import numpy as np
from astropy.table import Table


def test_example_acf():
    # DEFINE PARAMETERS  ==========================================================
    par = gun.packpars(kind="acf")
    par.doboot = False  # Do bootstrap error estimates
    par.dsept = 0.10
    par.nsept = 33
    par.septmin = 0.01
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
    cnt = gun.acf(gals, rans, par, nthreads=nt, plot=False, write=False)
    # ==============================================================================

    r_binedges_acf = np.load("tests/data/r_binedges_acf.npy")

    c_bin_acf = np.load("tests/data/c_binedges_acf.npy")
    l_binedges_acf = np.load("tests/data/l_binedges_acf.npy")
    rr_acf = np.load("tests/data/rr_acf.npy")
    dd_acf = np.load("tests/data/dd_acf.npy")
    w_acf_nat = np.load("tests/data/w_acf_nat.npy")

    assert len(c_bin_acf) == par.nsept, f"Expected {num_bins + 1} bin edges, but got {len(bin_edges)}"
    np.testing.assert_almost_equal(
        c_bin_acf, cnt.thm, decimal=9, err_msg="Bin centers do not match expected values"
    )
    np.testing.assert_almost_equal(
        r_binedges_acf, cnt.thr, decimal=9, err_msg="Bin right borders do not match expected values"
    )
    np.testing.assert_almost_equal(
        l_binedges_acf, cnt.thl, decimal=9, err_msg="Bin left borders do not match expected values"
    )
    assert isinstance(cnt.wth, np.ndarray), f"Expected hist to be a numpy array, but got {type(hist)}"
    assert isinstance(cnt.rr, np.ndarray), f"Expected hist to be a numpy array, but got {type(hist)}"
    assert isinstance(cnt.dd, np.ndarray), f"Expected hist to be a numpy array, but got {type(hist)}"
    assert np.issubdtype(cnt.wth.dtype, np.float), "cnt.wth dtype is not float"
    assert np.issubdtype(cnt.rr.dtype, np.float), "cnt.rr dtype is not float"
    assert np.issubdtype(cnt.dd.dtype, np.float), "cnt.dd dtype is not float"
    np.testing.assert_almost_equal(
        w_acf_nat, cnt.wth, decimal=9, err_msg="Correlation function is not correct"
    )
    np.testing.assert_almost_equal(
        rr_acf, cnt.rr, decimal=9, err_msg="Random-Random histogram is not correct"
    )
    np.testing.assert_almost_equal(
        dd_acf, cnt.dd, decimal=9, err_msg="Object-Object histogram is not correct"
    )
