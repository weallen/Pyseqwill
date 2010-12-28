from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
import os.path

numpy_path = os.path.join(numpy.__path__[0], 'core', 'include')

setup(
    name = "Pyseqwill",
    cmdclass = {'build_ext' : build_ext},
    ext_modules = [Extension("cHMM", ["cHMM.pyx"], include_dirs=[numpy.get_include()])]
)
