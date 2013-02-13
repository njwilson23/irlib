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


irview
======

*irview* is a graphical utility intended for viewing radar lines. Filters can
be applied or removed on the fly, and features within the radargram can be
manually digitized using a mouse.

[incomplete]


icepick and icerate
===================

*icepick* is a tool for precisely picking radar reflections (e.g. for bed
sounding). *icerate* is a similar program for rating the quality of picks
before surface interpolation. Both have an interface that is similar to that of
*irview*.

[incomplete]


irlib API
=========

The following sections describe the classes of the ``irlib`` API. These can
then be used directly from a Python script or terminal to view and manipulate
radar data.

Surveys
-------

.. automodule:: irlib.survey
    :members:

Gathers
-------

.. automodule:: irlib.gather
    :members:

Metadata
--------

.. automodule:: irlib.recordlist
    :members:


Management of picking and rating files
--------------------------------------

.. automodule:: irlib.filehandler
    :members:


itools
------

.. automodule:: irlib.itools
    :members:





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

