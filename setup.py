from distutils.core import setup
import numpy

from Cython.Build import cythonize

setup(
    name='test',
    ext_modules=cythonize("test.pyx"),
    include_dirs=[numpy.get_include()]
)
