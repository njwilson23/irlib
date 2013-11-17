from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext
from Cython.Build import cythonize

# Run with `python setup.py build_ext --inplace`

ext_modules = [Extension("irlib.agc", ["irlib/agc.pyx"])]
ext_modules = cythonize(ext_modules)

setup(
    name = "irlib",
    author = "Nat Wilson",
    packages = ["irlib"],
    ext_modules = ext_modules,
    #cmdclass = {'build_ext': build_ext},
    scripts = ["icepick2.py", "icerate.py", "h5_add_utm.py", "h5_dumpmeta.py", \
               "h5_replace_gps.py", "h5_generate_caches.py", "h5_consolidate.py", \
               "h52mat.py"],
)
