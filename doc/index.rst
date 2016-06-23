.. radar_tools documentation master file, created by
   sphinx-quickstart on Wed Feb 13 12:02:03 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for IcePick and irlib (radar_tools)
=================================================

.. toctree::
    :depth: 3

Introduction
============

.. image:: front_image.png


*radar_tools* is a Python package and set of applications that I wrote in order
to view and analyze radar data. I have made the package open-source
(http://njwilson23.github.com/radar_tools/), and written this documentation to
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

Dependencies
------------

*radar_tools* is built upon a number of standard tools from the scientific
Python ecosystem. The following are *required*:

- Python_ : Already installed for Linux/Mac OS X users

- Numpy_ : Basic array type, analogous to a matrix in MATLAB, except better

- Scipy_ : Wrappers for scientific libraries used for efficient filtering

- h5py_ : interface for HDF datasets

- matplotlib_ : Plotting library required for GUI tools

- pyproj_ : Wrapper for *proj.4* geographical projection library

Finally, these are *nice to have*:

- Cython_ : Python compiler for improving performance

- IPython_ : Interactive computing environment similar to MATLAB or Mathematica

Consider using a package manager (e.g. APT, rpm, pacman, or Homebrew).


Installation
------------

The latest version is on Github_. After downloading either directly or using the
command

::

    >> git clone git@github.com:njwilson23/irlib.git

Installation is best done with ``pip``, the Python package manager.

::

    >> cd irlib/    # [or whereever it's downloaded to]
    >> pip install .

Assuming that dependencies are available (see above), this will take care of
installing ``radar_tools`` properly. 

To use the *pywavelet* wavelet transform algorithms, navigate to
``irlib/external`` and follow the directions in the ``README`` file, being sure
to move the created file ``pywavelet.so`` to some place from which it can be
imported.

Installing manually
~~~~~~~~~~~~~~~~~~~

Alternatively, *irlib* can be build in place without ``pip`` by doing

::

    >> python setup.py build_ext --inplace

Path
~~~~

For convenience, programs that make up *radar\_tools* should be on the execution
``PATH``. If ``pip`` was used, this should be taken care of. Otherwise, on Linux
and Mac OS X, one can add the following line to the ``.bashrc``:

::

    export PATH=$PATH:~/python/irlib

On Windows, one should be able to modify the *Path* variable by right clicking
on **My Computer** and going to *Properties -> Advanced System Settings ->
Environment Variables*.

.. include:: doc_utilities.rst

.. include:: doc_icepick.rst

.. include:: doc_custom_commands.rst

.. include:: doc_tutorial.rst

.. include:: doc_irlib.rst

.. include:: doc_references.rst

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

