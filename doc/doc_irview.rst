irview
======

*irview* is a graphical utility intended for viewing radar lines. Filters can
be applied or removed on the fly, and features within the radargram can be
manually digitized using a mouse.

Usage
-----

*irview* is invoked from the command line using the syntax

::

    irview.py -f HDFNAME [-L LINENO]

where ``HDFNAME`` is the name of the radar survey file and ``LINENO`` is an
optional line number (default ``0``).

Upon launch, the terminal window is converted into an *irview* console, and a
figure window showing the radargram (B-scan) opens. Clicking on the radargram
prints information about the point clicked in the console. For example:

::

        FID: 0005004700000000
        x: 26       y:188           t: 1503.61 ns

means that the point clicked is the 26th displayed trace at the 188th sample,
which is approximately 1503.61 ns from the top of the radargram. Each trace is
assigned a *FID*, which is a 16 character string indicating:

.. table:: FID meaning

    =========== ============================
    Chars       Description
    =========== ============================
     1-4        Line number
     5-8        Trace number
     9-12       Datacapture (channel) number
     13-16      Echogram number (unused)
    =========== ============================

The fact that the trace number (``47``) according to the *FID* is not the same
as the display number of the trace (``26``) arises because some traces have
been filtered out of the dataset prior to display (perhaps because they didn't
trigger properly, or very close to other traces). The *FID* is intended to be a
static identifier permanently bound to the original data, and is not recomputed
when a dataset is modified.

Commands can be typed in the console. Typing ``help`` gives an idea of the
options. Important examples include:

- ``info`` provides more detailed information about the currently opened line. 

- ``open LINENO`` allows different lines to be viewed without restarting *irview*

- ``exit`` closes *irview*.


Caching
-------

The performance of *irview* (and any other program in *radar_tools* that loads
radar lines) can be enhanced substantially by pre-caching of the radar lines.

[**incomplete - should probably be a separate chapter**]


Filters
-------

In addition to viewing raw data, *irview* provides access to pre-defined
filters on the fly using the ``f`` command (mnemonic "filter"). For example,
typing::

    f gc

applies a linear gain control operator to each trace, increasing the amplitude
of later-arriving events. This modifies the data in memory, and **it does not
alter the data in the HDF file**. In order to determine what filters have been
applied to a dataset, typing ``f`` alone lists them, along with any parameters
required to reproduce the presently displayed data. The data can be reset to
the version originally loaded by typing ``nf`` (mnemonic "no filter").

The presently defined filters include common operations such as time-domain and
frequency-domain filtering, dewow, linear and automatic gain control, F-K
migration, and instrument ringing removal. Defining custom filters requires
modification of the file ``irlib/filter_defs.py``, and can be done using the
*irlib* API (discussed elsewhere).


Digitizing
----------

A tangential application of *irview* is hand-digitizing of radar features. This
is less precise that trace-by-trace picking (see *icepick*), but more
appropriate for selecting volumetric features, or features for which the
individual traces are messy (geological applications?).

To begin digitizing a feature, type::

    dnew [FEATURENAME]
    
where ``FEATURENAME`` is an optional descriptive word or comment. The radargram
can then be clicked with the left mouse button to delineate shapes by vertex.
Pressing the middle mouse button with "undo" the last vertex created. Pressing
the right mouse button will create a final vertex and then end the feature.

*Alternatively, pressing "N" (Shift+n) while the figure window is focused can
be used to start a new feature (with no comment) and "E" (Shift+e) will end the
feature.*

Once all desired features have been digitized, typing::

    dexport

saves the vertices to a text file. Presently, there is no way to load
previously-created features across sessions.

Additional commands:

- ``dls`` lists previously-created features

- ``drm NUMBER`` deletes the feature identified by ``NUMBER``




