# -*- coding: utf-8 -*-
"""
Example code to compare two pre-made correlations (red and blue galaxies)
against a control correlation (all galaxies) and fitt two powerlaws

Note :
   
Make sure to first unzip ./samples_data/red_blue_all.zip, which contains
the 3 .cnt files
"""
import matplotlib.pyplot as plt
import gundam as gun
import os, sys

# DEFINE DATA
fred  = './examples_data/red.cnt'
fblue = './examples_data/blue.cnt'
fall  = './examples_data/all.cnt' 

if (not(os.path.isfile(fred)) or not(os.path.isfile(fblue)) or not(os.path.isfile(fall))):
    sys.exit('Make sure to unzip first "./samples_data/red_blue_all.zip"')
# MAKE COMPARISON CF PLOT (no need to read .cnt explicitly)
f, ax1, ax2 = gun.comparecf([fred, fblue, fall], [fall], plotratio=True, fill=True)

# NOW READ CORRELATIONS
red  = gun.readcounts(fred)
blue = gun.readcounts(fblue)

# FIT POWER LAW TO RED AND BLUE GALAXIES BETWEEN 0.1-10 Mpc
plt.sca(ax1)
gun.fitpowerlaw(red.rpm, red.wrp, red.wrperr, fitrange=[0.1,10.], plot=True)
gun.fitpowerlaw(blue.rpm, blue.wrp, blue.wrperr, fitrange=[0.1,10.], plot=True)

plt.show()

