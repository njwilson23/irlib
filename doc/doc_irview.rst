irview
======

*irview* is a graphical utility intended for viewing radar lines. Filters can
be applied or removed on the fly, and features within the radargram can be
manually digitized using a mouse. All GUI programs require *matplotlib* (see
Dependencies_).

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
This can be done using the API (``Gather.Dump()``), or by running the
commandline utility ``h5_generate_caches`` (discussed previously). Caches
are Python "pickles" (serialized data), and contain a snapshot of the radar
data, as well as a reference to ``irlib``. Substantial changes to ``irlib``
sometimes require re-caching.


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

**Non-comprehensive list of filters**

=============== ==============================================================
Command         Description
=============== ==============================================================
*Gain Control*
``gc``          Applies a linear gain enhancement
``gc2``         Applies a quadratic gain enhancement
``agc``         Applies an automatic (nonlinear) gain enhancement
*Covolutions*
``dewow``       "Dewowing" filter to remove instrument drift
``lowpass``     Performs a frequency-domain lowpass filter with a cutoff
                frequency of 25 MHz
``highpass``    Performs a frequency-domain lowpass filter with a cutoff
                frequency of 25 MHz
``lowpass_ma``  Performs a time-domain lowpass filter (moving average)
``highpass_ma`` Performs a time-domain highpass filter (inverted moving
                average)
*Recursive*
``iir30low``    Chebyschev lowpass filter with cutoff at 30 MHz
``iir25high``   Chebyschev highpass filter with cutoff at 25 MHz
*Migration*
``fkmig``       Stolt (F-K) migration
*Misc*
``abs``         Displays the absolute value of the data
``wiener``      Wiener statistical noise filter
``ringing``     Horizontal ringing filter based on singular value decomposition
``project``     Project radar line to straight segments with equal trace
                spacing
=============== ==============================================================

For a comprehensive list, see ``irlib/filter_defs.py:ApplyFilter()``


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

    dsave

saves the vertices to a text file. The saved file is Tab-delimited with blank
lines between features.

=========== ===================================
Column      Description
=========== ===================================
1           Trace at vertex location
2           Longtitude
3           Latitude
4           Time (ns) from the top of the trace
=========== ===================================

Presently, comments are not saved in the file, and there is no way to load
previously-created features across sessions.

Additional commands:

- ``dls`` lists previously-created features

- ``drm NUMBER`` deletes the feature identified by ``NUMBER``




