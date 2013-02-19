icepick and icerate
===================

*icepick* is a tool for precisely picking radar reflections (e.g. for bed
sounding). *icerate* is a similar program for rating the quality of picks
before surface interpolation. Both have an interface that is similar to that of
*irview*, so this will assume that the previous chapter has been read and focus
on the differences from *irview*.

icepick
-------

The syntax for opening *icepick* and *icerate* is similar to *irview*::

    icepick.py -f HDFNAME [-L LINENO]

The window that opens differs slightly. In *icepick*, the radargram is
complemented by a large A-scan view that shows individual radar traces. The
location along the radargram that the traces are from is shown by the vertical
yellow lines in the radargram pane. Assuming there are more than eight traces
(normally the case), the display can be panned across the radargram with the
**h** and **l** (ell) keys (which should be second-nature to users of *Vi*).

The process of picking is fairly simple. In the lower panel of the icepick
window (where the individual traces are shown), click the mouse on the part of
the trace representing a reflection to be timed. Right-clicking removes the
pick if you've made a mistake. Fine adjustments can be made by pressing the
**j** (down) and **k** (up) keys. Whenever the side-scrolling keys are pressed
(**h** and **l**), a line representing the picks is drawn on the radargram.
Presumably, the bed should be picked on every trace where it can be identified.

Filters can be applied as in *irview*. An additional tool to assist with
picking is the ``autobed [min, max]`` command, which tries to identify the bed
reflection wavelet (between samples ``min`` and ``max``; see the right vertical
axis) based on a set of simple heuristics. This works well when the bed
reflection is obvious and the radargram is clean, but will typically require
manual clean-up.

The "picking mode" can be switched between bed-picking and airwave-picking by
typing ``mode x`` where ``x`` is either ``bed`` or ``dc``. If ``x`` is not
specified, then the current mode will be printed to the console.  There is an
``autodc`` command that is analogous to ``autobed``, except that it tends to be
much more successful.

Once the picks are satisfactory, type ``save`` to save the timing data to the
folder ``picking/``.

icerate
-------

Once picks have been made, they can be rated for quality. *icerate* shows the
picks made in *icepick* (either in order or randomly), and allows a numeric
quality rating to be applied (individually or in small groups) by typing
**1-5**.

The interface of *icerate* should be familiar after *irview* and *icepick*.
Like *icepick*, there is an A-scan view, but it is smaller and out of the way.
When ``save`` is typed, the results will go into ``rating/``.

