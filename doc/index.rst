.. radar_tools documentation master file, created by
   sphinx-quickstart on Wed Feb 13 12:02:03 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for radar_tools
=============================

.. toctree::
    :depth: 3

Introduction
============

*radar_tools* is a Python package and set of applications that I wrote in order
to view and analyze radar data. I have made the package open-source
(https://github.com/njwilson23/radar_tools), and written this documentation to
make it easier for others to use.

There are two sides to *radar_tools*. The first is the set of command-line and
graphical utilities that perform common operations on ice-penetrating radar
data. This includes concatenating datasets, projecting spatial coordinates,
extracting metadata, viewing and filtering radar lines, and picking reflection
wavelets. A tutorial following this half of the manual steps through a typical
workflow for processing radar data. For the Python-doubtful, a script is
included in this section that converts raw HDF datasets into matfiles for
processing in MATLAB.

The second part of *radar_tools* is the ``irlib`` API, which is used by
creating custom Python scripts. ``irlib`` serves the role of abstracting radar
data into easily-used datastructures with built-in analysis functionality.
Accessing the API directly is useful for experimental data exploration, filter
construction, individualized plotting, and interfacing *radar_tools*
capabilities with external Python scripts. The commandline and GUI utilities in
*radar_tools* are themselves built on ``irlib``.

Getting the code
----------------

Although a version of *radar_tools* will be present on the SFU Lab NAS, I've
been updating the source code fairly frequently, and changes nearly always make
it better. I recommend pulling the latest version from Github_. The ideal way
to do this is to use git_. This way, new changes can be fetched with a single
command rather than re-downloading and unpacking a compressed archive, and
fixes that you happen to make can be easily shared upstream. Alternatively
one-off zipped snapshot can be downloaded and installed.

Dependencies
------------

*radar_tools* is built upon a number of standard tools from the scientific
Python ecosystem. The following are *required*:

- Python_ : Already installed for Linux/Mac OS X users

- Numpy_ : Basic array type, analogous to a matrix in MATLAB, except better

- Scipy_ : Wrappers for scientific libraries used for efficient filtering

- h5py_ : interface for HDF datasets

The following are *required for individual tools to work*:

- matplotlib_ : Plotting library required for GUI tools

- pyproj_ : Wrapper for *proj.4* geographical projection library

Finally, these are *nice to have*:

- Cython_ : Python compiler for improving performance

- IPython_ : Interactive computing environment similar to MATLAB or Mathematica

Consider using a package manager (e.g. APT, rpm, pacman, or Homebrew).

Installation and set-up
-----------------------

Basic installation is as simple as ensuring the required Dependencies_ are met
and either ``git clone``-ing or unzipping *radar_tools*. In order to use the
compiled AGC function, type

::

    python setup.py build_ext --inplace

To use the *pywavelet* wavelet transform algorithms, navigate to
``irlib/external`` and follow the directions in the ``README`` file, being sure
to move the created file ``pywavelet.so`` to some place from which it can be
imported.

Path
~~~~

For convenience, programs that make up *radar\_tools* should be on the
execution ``PATH``. On Linux and Mac OS X, one can add the following line to my
``.bashrc``:

::

    export PATH=$PATH:~/python/irlib

On Windows, one should be able to modify the *Path* variable by right clicking
on **My Computer** and going to *Properties -> Advanced System Settings ->
Environment Variables*.



.. include:: doc_utilities.rst

.. include:: doc_irview.rst

.. include:: doc_icepick.rst

.. include:: doc_tutorial.rst

.. include:: doc_irlib.rst

.. include:: doc_references.rst

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

