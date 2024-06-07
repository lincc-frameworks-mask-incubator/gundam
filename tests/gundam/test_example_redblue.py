"""
Example code to compare two pre-made correlations (red and blue galaxies)
against a control correlation (all galaxies) and fitt two powerlaws
"""
from gundam import gundam as gun

# DEFINE DATA
fred = "tests/data/red_blue_all/red.cnt"
fblue = "tests/data/red_blue_all/blue.cnt"
fall = "tests/data/red_blue_all/all.cnt"


def test_example_redblue():
    # MAKE COMPARISON CF PLOT (no need to read .cnt explicitly)
    f, ax1, ax2 = gun.comparecf([fred, fblue, fall], [fall], plotratio=True, fill=True)

    # NOW READ CORRELATIONS
    red = gun.readcounts(fred)
    blue = gun.readcounts(fblue)

    # FIT POWER LAW TO RED AND BLUE GALAXIES BETWEEN 0.1-10 Mpc
    # plt.sca(ax1)
    gun.fitpowerlaw(red.rpm, red.wrp, red.wrperr, fitrange=[0.1, 10.0], plot=True)
    gun.fitpowerlaw(blue.rpm, blue.wrp, blue.wrperr, fitrange=[0.1, 10.0], plot=True)
    # plt.show()
