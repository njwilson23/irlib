irlib API
=========

The following sections describe the classes of the ``irlib`` API. These can
then be used directly from a Python script or terminal to view and manipulate
radar data.

Surveys (Collections of lines in HDF Files)
-------------------------------------------

.. automodule:: irlib.survey
    :members:

Gathers (Common-offset, Common-midpoint, etc.)
----------------------------------------------

.. automodule:: irlib.gather
    :members:

Metadata
--------

.. automodule:: irlib.recordlist
    :members:


Picking and rating file management
----------------------------------

.. automodule:: irlib.filehandler
    :members:


itools
------

.. automodule:: irlib.itools
    :members:

Application building
--------------------

The following modules contain the building blocks the graphical applications
(such as IcePick2), as well as user-defined filters.

Console
~~~~~~~

.. automodule:: irlib.app.console
    :members:

Windows
~~~~~~~

.. automodule:: irlib.app.components
    :members:

Command parsing
~~~~~~~~~~~~~~~

.. automodule:: irlib.app.command_parser
    :members:

Filter framework
~~~~~~~~~~~~~~~~

.. automodule:: irlib.app.filters
    :members:



