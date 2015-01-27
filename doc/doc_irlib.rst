irlib API
=========

The following sections describe the classes of the ``irlib`` API. These can
then be used directly from a Python script or terminal to view and manipulate
radar data.

Surveys (Collections of lines in HDF Files)
-------------------------------------------

.. automodule:: irlib.survey

.. autoclass:: irlib.survey.Survey
    :members:

.. autoclass:: irlib.survey.EmptyLineError

Gathers (Common-offset, Common-midpoint, etc.)
----------------------------------------------

.. automodule:: irlib.gather

.. autoclass:: irlib.gather.Gather
    :members:

.. autoclass:: irlib.gather.CommonOffsetGather
    :members:

.. autoclass:: irlib.gather.CommonMidpointGather
    :members:

.. autoclass:: irlib.gather.PickableGather
    :members:

.. autoclass:: irlib.gather.LineGatherError

Metadata
--------

.. automodule:: irlib.recordlist

.. autoclass:: irlib.recordlist.RecordList
    :members:

.. autoclass:: irlib.recordlist.ParseError

Picking and rating file management
----------------------------------

.. automodule:: irlib.filehandler

.. autoclass:: irlib.filehandler.FileHandler
    :members:

.. autoclass:: irlib.filehandler.FileHandlerError


itools
------

.. automodule:: irlib.itools
    :members:

Application building
--------------------

The following modules contain the building blocks the graphical applications
(such as IcePick2), as well as user-defined filters. They would be used to
develop or extent the existing graphical tools.

Console
~~~~~~~

.. automodule:: irlib.app.console

.. autoclass:: irlib.app.console.Console
    :members:

Windows
~~~~~~~

.. automodule:: irlib.app.components

.. autoclass:: irlib.app.components.AppWindow
    :members:

.. autoclass:: irlib.app.components.Radargram
    :members:

.. autoclass:: irlib.app.components.PickWindow
    :members:

.. autoclass:: irlib.app.components.MapWindow
    :members:

Command parsing
~~~~~~~~~~~~~~~

.. automodule:: irlib.app.command_parser
    :members:

Filter framework
~~~~~~~~~~~~~~~~

.. automodule:: irlib.app.filters
    :members:



