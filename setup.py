from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Build import cythonize
    ext_modules = [Extension("irlib.agc", ["irlib/agc.pyx"])]
    ext_modules = cythonize(ext_modules)
except ImportError:
    ext_modules = []

setup(
    name = "irlib",
    author = "Nat Wilson",
    version = "0.4",
    packages = ["irlib"],
    ext_modules = ext_modules,
    scripts = ["icepick2.py", "icerate.py", "h5_add_utm.py", "h5_dumpmeta.py", \
               "h5_replace_gps.py", "h5_generate_caches.py", "h5_consolidate.py", \
               "h52mat.py"],
)
