#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# import os
# from setuptools import setup
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

# f2py -c cflibfor.pyf cosmolib.f90 cflibfor.f90


def configuration(parent_package="", top_path=None):
    config = Configuration(None, parent_package, top_path)
    libraries = [
        "gomp",
        # 'blas',
    ]
    config.add_extension(
        "cflibfor",
        ["cflibfor.pyf", "cflibfor.f90"],
        libraries=libraries,
        # f2py_options = ['--debug-capi'], #'--debug-capi'
        # define_macros = [('F2PY_REPORT_ON_ARRAY_COPY','1')],
        # this is the flag gfortran needs to process OpenMP directives
        # -fopt-info-vec-missed  -funsafe-math-optimizations #
        # -ftree-vectorizer-verbose=6  -fopt-info-all...
        # -ffast-math -march=corei7-avx
        extra_compile_args=["-fopenmp"],
        extra_f90_compile_args=[
            "-march=native",
            "-ftree-vectorize",
            "-funroll-loops",
        ],
        # extra_link_args = ['-dynamic'],
    )
    return config


if __name__ == "__main__":
    setup(
        name="gundam",
        version="0.81",
        configuration=configuration,
        packages=["gundam"],
    )
