################################################################################
# Example of building Gundam with distutils and compile automatically with f2py
#
# Just type "make" to build the .so library
#
# If this does not work, try compiling directly with f2py :
#
#  -> f2py -c --opt='-O3 -march=corei7-avx -ftree-vectorize' cflibfor.pyf cflibfor.f90
#
################################################################################

.PHONY: all clean test

all:
#	python setup.py build_ext --inplace --fcompiler=gnu95
	env F90FLAGS=-fopenmp python setup.py build_ext --inplace --fcompiler=gnu95

test: all
	@echo Running gundam test...yet TODO
	python script.py

clean:
	rm -rf build
	rm -f cflibfor-f2pywrappers2.f90
	rm -f cflibformodule.c
	rm -f cflibfor*.so
