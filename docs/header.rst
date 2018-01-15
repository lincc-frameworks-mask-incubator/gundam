============================================================
GUNDAM :  A Toolkit For Fast Two-Point Correlation Functions
============================================================

Gundam is a package to **count pairs** and **calculate spatial two-point correlation 
functions** of large galaxy samples. Its main features are :

Speed
    By calling Fortran routines that implement efficient skip-list/linked-list
    algorithms, Gundam can be extremely fast
  
Parallel
    Can automatically run in parallel to use all cores available. It employs the
    `OpenMP <http://www.openmp.org>`_ framework to make use of multi-core CPUs

User-friendly
    By carefully wrapping Fortran code in a suitable Python framework, Gundam is 
    very easy to use. A typical run consists of just 3 lines of code : (1) read 
    data, (2) define parameters, (3) get counts

Error estimates, user-defined weights, fiber corrections
    Gundam can estimate boostrap errors, weight pair counts, and even correct 
    counts for fiber collisions

Plotting functions
    Gundam can produce nice, paper ready plots for 1D and 2D correlations,
    complete with ratios, labels and even power-law fits
    
Extensible
    Desgined in 3 layers of main, auxiliary and wrapper routines, it is quite
    easy to extend functionality by novice as well as seasoned users

    
Pair counts and correlation functions can be saved in ASCII files, as well as 
in a dictionary-like object that holds all calculations, input parameters 
and log messages. Share this object with your collaborators instead of just
the final plot.

Though intented primarly for redshift surveys, it can also be adapted for 
simulation data and ultimately for any set of points in space.


@author: Emilio Donoso <edonoso@conicet.gov.ar>
