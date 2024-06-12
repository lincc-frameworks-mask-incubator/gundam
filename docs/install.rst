.. _step_by_step_install:

************************
Installation
************************

To install Gundam, you have two choices: (1) use pip, or (2) build from scratch.
Method (1) is useful if you simply want to use gundam. Method (2) on the other
hand allows easy access to modify or extend the Fortran counting routines.

Option 1: Building from pip
==============================

To install via pip execute:

.. code:: python

      pip install gundam

Option 2: Building from source
=========================================

If you want to contribute to the package, you need to clone the Gundam repository and
install with the "dev" optional dependencies.

.. code::

      git clone https://github.com/lincc-frameworks-mask-incubator/gundam.git
      cd gundam
      pip install -e .'[dev]'
      pre-commit install

By default, this will compile and build the library in-place. Feel free to modify
the `CMakeLists` file to suit your needs.

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

