.. icepick and irlib (radar_tools) documentation master file, created by
   sphinx-quickstart on Sat Sep 18 06:50:40 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for IcePick and irlib (radar_tools)
=================================================

.. toctree::
   :maxdepth: 3
   :caption: Contents:
   
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

Irlib was first developed prior to 2012 back in Python2 days.  Since then Python
2 is no longer supported so everything from version 0.5 on works with Python 3.  
Some of the functions may be backward compatible but this has not been fully 
tested.  If you need Python 2 functionality, please use and older version of 
irlib.  Users may also wish to check out another open source suite of radar 
tools.  See: https://impdar.readthedocs.io/en/latest/ 

The following documentation is as complete as possible.  For the **beginner user** 
it starts with how to install the software and then progresses through some of the 
data handling utilities and then reviews how to pick radar returns.  Following 
this, an **intermediate user** might wish to follow the in-depth tutorial.  
At the end of the manual there are topics for **advanced users and developpers**. 

**Beginner User**
Introduction
Installation
Command-line Utilities
IcePick2
icerate

**Intermediate User**
Tutorial

**Advanced Users and Developpers**
Changes in IEI h5 file format
Adding customized filters to the irlib GUI apps
Documentation
irlib API

.. introduction is above... next installation

.. include:: doc_installation.rst

.. Commandline utilities section
.. include:: doc_utilities.rst

.. Icepick 
.. include:: doc_icepick.rst

.. Tutorial is more in-depth and 
.. include:: doc_tutorial.rst

.. include:: doc_h5formats.rst

.. include:: doc_custom_commands.rst

.. Sphinx documentation 
.. include:: doc_documentation.rst

.. irlib API section
.. include:: doc_irlib.rst  


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

