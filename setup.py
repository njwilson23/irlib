from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# Run with `python setup.py build_ext --inplace`

ext_modules = [Extension("irlib.agc_cy", ["irlib/agc.pyx"])]

setup(
  name = 'Fast AGC function',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
