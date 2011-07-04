#!/usr/bin/env python
#
# Distutils build file
#

from sys import argv
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

##############################################################################

if __name__ == "__main__":
    setup(
        name='sparse',
        version='0.1',
        description='',
        author='',
        author_email='',
        license='TBD',
        cmdclass = {'build_ext': build_ext},
        ext_modules = [
            Extension(
                "sparse.wrapper",
                ["sparse/wrapper.pyx"],
                language="c++",
            )
        ],
        packages=['sparse']
    )

