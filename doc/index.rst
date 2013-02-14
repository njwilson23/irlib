.. radar_tools documentation master file, created by
   sphinx-quickstart on Wed Feb 13 12:02:03 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for radar_tools
=============================

.. contents::
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
wavelets.

The second part is the ``irlib`` API, upon which the first is built. Accessing
the API directly is useful for experimental data exploration, filter
construction, and interfacing *radar_tools* capabilities with external Python
scripts.



.. irview
.. include:: doc_irview.rst


.. icepick and icerate
.. include:: doc_icepick.rst


.. irlib API
.. include:: doc_irlib.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

