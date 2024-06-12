.. _step_by_step_install:

************************
Installation
************************

To install Gundam, you have two choices: (1) build from scratch, or (2) use pip. 
I recommend method (1), since it will allow easy access to modify or extend the 
Fortran counting routines. In any case, make sure to fulfill the required
:ref:`gundam_dependencies` before installation.

          
Option 1: Building from source
==============================

Just install the required dependencies, clone the Gundam repository and type `make`

.. code:: 
          
	  git clone https://github.com/lincc-frameworks-mask-incubator/gundam.git
	  cd gundam
	  make

By default this will compile and build the library in-place. Feel free to modify 
the `Makefile` to suit your needs. After compilation, you can optionally install
the library in your default global-site packages directory

.. code:: 

	  python setup.py install

or in your default user packages directory

.. code:: 

	  python setup.py install --user


Option 2: Using pip (not yet functional!)
=========================================

.. code:: python

          pip install gundam

          
.. _gundam_dependencies:
          
Dependencies
============
0. `Python <http://www.python.org/>`_: 3.5 or later 
1. `GCC compiler (Fortran & C) <https://gcc.gnu.org/>`_
2. `Numpy <http://www.numpy.org/>`_, `matplotlib <http://matplotlib.org/>`_, 
   `astropy <http://www.astropy.org/>`_, `scipy <https://www.scipy.org/>`_
3. `munch <https://pypi.python.org/pypi/munch>`_
4. `pymorton <https://github.com/trevorprater/pymorton/>`_ (only needed when playing with Morton ordering, instead of the default "pixel" sorting)

Any of the above can be installed with pip or conda. A few other codes for Morton
and Hilbert ordering are also included if you wish to experiment with alternative 
sorting methods, but are not really needed if you stick with the default 
sorting scheme).

Gundam has been tested under Anaconda 4.4.0, Python 3.6.1, Numpy 1.12.1 and gcc 7.1.1, running OpenSuse
in an Intel platform.

