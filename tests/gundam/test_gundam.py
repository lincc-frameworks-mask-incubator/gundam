import pytest
import gundam as gun
from astropy.table import Table
from pathlib import Path


@pytest.fixture(scope='module')
def input_data():
    local_dir = Path(__file__).resolve().parent
    galf = "DATA.fits"  # Galaxy sample
    ranf = "RAND.fits"  # Random sample
    crossf = "DR7-lrg.fits"
    galf_file_path = local_dir.parent / 'data' / galf
    ranf_file_path = local_dir.parent / 'data' / ranf
    cross_file_path = local_dir.parent / 'data' / crossf
    gals = Table.read(galf_file_path)
    if "wei" not in gals.colnames:
        gals["wei"] = 1.0  # If not present, set weights to 1
    rans = Table.read(ranf_file_path)
    if "wei" not in rans.colnames:
        rans["wei"] = 1.0  # If not present, set weights to 1
    cross = Table.read(cross_file_path)
    if "wei" not in cross.colnames:
        cross["wei"] = 1.0  # If not present, set weights to 1
    return gals, rans, cross


@pytest.mark.fast
def test_example_acf(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="acf")
    par.autogrid = False  # Automatic SK grid size
    par.mxh1 = 60  # SK size in dec
    par.mxh2 = 240  # SK size in ra
    par.mxh3 = 24  # SK size in z
    par.nsepp = 78  # Number of bins of projected separation rp
    par.seppmin = 0.01  # Minimum rp in Mpc/h
    par.dsepp = 0.5  # Bin size of rp (in log space)
    par.logsepp = 0  # Use linear bins instead of log bins
    par.nsepv = 60  # Number of bins of LOS separation pi
    par.doboot = False  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    cnt = gun.acf(gals, rans, par, nthreads=10, plot=True, write=False)


@pytest.mark.slow
def test_example_acf(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="acf")
    par.autogrid = False  # Automatic SK grid size
    par.mxh1 = 60  # SK size in dec
    par.mxh2 = 240  # SK size in ra
    par.mxh3 = 24  # SK size in z
    par.nsepp = 78  # Number of bins of projected separation rp
    par.seppmin = 0.01  # Minimum rp in Mpc/h
    par.dsepp = 0.5  # Bin size of rp (in log space)
    par.logsepp = 0  # Use linear bins instead of log bins
    par.nsepv = 60  # Number of bins of LOS separation pi
    par.doboot = False  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    cnt = gun.acf(gals, rans, par, nthreads=1, plot=True, write=False)


@pytest.mark.fast
def test_example_acf_ls(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="acf")
    par.autogrid = False  # Automatic SK grid size
    par.mxh1 = 60  # SK size in dec
    par.mxh2 = 240  # SK size in ra
    par.mxh3 = 24  # SK size in z
    par.nsepp = 78  # Number of bins of projected separation rp
    par.seppmin = 0.01  # Minimum rp in Mpc/h
    par.dsepp = 0.5  # Bin size of rp (in log space)
    par.logsepp = 0  # Use linear bins instead of log bins
    par.nsepv = 60  # Number of bins of LOS separation pi
    par.doboot = False  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "LS"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    cnt = gun.acf(gals, rans, par, nthreads=10, plot=True, write=False)


@pytest.mark.fast
def test_example_accf(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="accf")
    par.autogrid = False  # Automatic SK grid size
    par.mxh1 = 60  # SK size in dec
    par.mxh2 = 240  # SK size in ra
    par.mxh3 = 24  # SK size in z
    par.nsepp = 78  # Number of bins of projected separation rp
    par.seppmin = 0.01  # Minimum rp in Mpc/h
    par.dsepp = 0.5  # Bin size of rp (in log space)
    par.logsepp = 0  # Use linear bins instead of log bins
    par.nsepv = 60  # Number of bins of LOS separation pi
    par.doboot = False  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "NAT"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    par.cra2 = 'ra'
    par.cdec2 = 'dec'
    # CALCULATE THE CORRELATION
    cnt = gun.accf(gals, rans, cross, par, nthreads=10, write=False, plot=True)


@pytest.mark.fast
def test_example_accf_ls(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="accf")
    par.autogrid = False  # Automatic SK grid size
    par.mxh1 = 60  # SK size in dec
    par.mxh2 = 240  # SK size in ra
    par.mxh3 = 24  # SK size in z
    par.nsepp = 78  # Number of bins of projected separation rp
    par.seppmin = 0.01  # Minimum rp in Mpc/h
    par.dsepp = 0.5  # Bin size of rp (in log space)
    par.logsepp = 0  # Use linear bins instead of log bins
    par.nsepv = 60  # Number of bins of LOS separation pi
    par.doboot = False  # Do bootstrap error estimates
    par.omegam = 0.25  # Omega matter
    par.omegal = 0.75  # Omega lambda
    par.h0 = 100  # Hubble constant [km/s/Mpc]
    par.calcdist = True  # Calculate comov. dist.
    par.estimator = "LS"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    par.cra2 = 'ra'
    par.cdec2 = 'dec'
    # CALCULATE THE CORRELATION
    cnt = gun.accf(gals, rans, cross, par, nthreads=10, write=False, plot=True)


@pytest.mark.fast
def test_example_pcf(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="pcf")
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
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cred = 'z_spec'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    par.cred1 = 'z'
    # CALCULATE THE CORRELATION
    cnt = gun.pcf(gals, rans, par, nthreads=10, write=False, plot=True)


@pytest.mark.fast
def test_example_pcf_ls(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="pcf")
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
    par.estimator = "LS"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cred = 'z_spec'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    par.cred1 = 'z'
    # CALCULATE THE CORRELATION
    cnt = gun.pcf(gals, rans, par, nthreads=10, write=False, plot=True)


@pytest.mark.fast
def test_example_pccf(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="pccf")
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
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cred = 'z_spec'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    par.cred1 = 'z'
    # CALCULATE THE CORRELATION
    cnt = gun.pccf(gals, rans, cross, par, nthreads=10, write=False, plot=True)


@pytest.mark.fast
def test_example_pccf_ls(input_data):
    gals, rans, cross = input_data
    par = gun.packpars(kind="pccf")
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
    par.estimator = "LS"  # Choose Landy-Szalay estimator for the PCF
    par.cra = 'ra_coord'
    par.cdec = 'dec_coord'
    par.cred = 'z_spec'
    par.cra1 = 'ra'
    par.cdec1 = 'dec'
    par.cred1 = 'z'
    # CALCULATE THE CORRELATION
    cnt = gun.pccf(gals, rans, cross, par, nthreads=10, write=False, plot=True)


@pytest.mark.fast
def test_example_redblue():
    gundam_dir = Path(__file__).resolve().parent
    red_file = "red.cnt"
    blue_file = "blue.cnt"
    all_file = "all.cnt"
    red_dir = gundam_dir.parent  / 'data/red_blue_all' / red_file
    blue_dir = gundam_dir.parent  / 'data/red_blue_all' / blue_file
    all_dir = gundam_dir.parent  / 'data/red_blue_all' / all_file
    # MAKE COMPARISON CF PLOT (no need to read .cnt explicitly)
    f, ax1, ax2 = gun.comparecf([red_dir, blue_dir, all_dir], [all_dir], plotratio=True, fill=True)
    # NOW READ CORRELATIONS
    red = gun.readcounts(red_dir)
    blue = gun.readcounts(blue_dir)
    # FIT POWER LAW TO RED AND BLUE GALAXIES BETWEEN 0.1-10 Mpc
    # plt.sca(ax1)
    gun.fitpowerlaw(red.rpm, red.wrp, red.wrperr, fitrange=[0.1, 10.0], plot=True)
    gun.fitpowerlaw(blue.rpm, blue.wrp, blue.wrperr, fitrange=[0.1, 10.0], plot=True)
    # plt.show()
