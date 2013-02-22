Tutorial
========

I am going to use the radar data collected in Spring 2012 to create a bed map
for Glacier 3 ("East Glacier"). I'll document everything I do here, so that
this can serve as a step-by-step tutorial on on using the radar interpretation
tools. I'm working from a Linux computer, so adapt terminal commands as
necessary.

The radar codes are subject to change. I'll be using a current git revision as
of February 13, 2013 (``commit 3dc4638de82867b58c40cd19727d3cfd980112f6``).
Most likely, it will be best to use the most recent version available.
Furthermore, I'll be using gstat_ for interpolation.

Set-up
------

File structure
~~~~~~~~~~~~~~

I've created a directory to work from. I named it ``gl3bedmap``, and I added a
subdirectory called data that contains a copy of the raw radar data. Create
subdirectories called ``cache``, ``picking``, ``rating`` and ``offsets`` too,
as they'll be needed later. The directory tree should now look like:

::

    gl3bedmap:
        data:
            gl3_radar_may14.h5
            gl3_radar_may16-17.h5
        cache:
        picking:
        rating:
        offsets:

Concatenating datasets
~~~~~~~~~~~~~~~~~~~~~~

The field data are contained in two HDF5 datasets. Since they both come from
the same glacier and the same field campaign, it would be nice to just deal
with one file. This is what the *h5\_consolidate.py* tool is for. From the
``data`` directory:

::

    h5_consolidate.py gl3_radar_may14 gl3_radar_may16 -o gl3_radar_2012.h5

which creates a single dataset with all of the data. The line numbers in all
but the first dataset are changed to avoid repeats. If I were to type the name
``h5_consolidate.py`` by itself, it would print out a reminder of how it works.

Viewing the metadata
~~~~~~~~~~~~~~~~~~~~

It's convenient to be able to read the data we're dealing with. There are a
number of ways to do this. HDF software comes with a number of utilities,
including *h5dump* and *hdfview*. I find these good for looking at raw data,
but the output can be a bit clunky to wade through. Therefore, I wrote
*h5\_dumpmeta*, which parses each radar sounding and writes a CSV file that can
be opened quickly in MATLAB (or as a spreadsheet).

::

    h5_dumpmeta.py data/gl3_radar_2012.h5 > radar.csv

Note that this is just the metadata, and doesn't actually print out the
receiver sled's digitizer readings.

FID Interpretation
~~~~~~~~~~~~~~~~~~

Every row in the metadata is assigned an identification (FID) that is unique
within the file. The FID is a sixteen character string.

+------------------+------------------+----------------------------+
| Characters       | Meaning          | Notes                      |
+==================+==================+============================+
| 1-4              | Line             | Collections of traces      |
+------------------+------------------+----------------------------+
| 5-8              | Trace            | Individual soundings       |
+------------------+------------------+----------------------------+
| 9-12             | Datacapture      | Used for channels in       |
|                  |                  | dualdar                    |
+------------------+------------------+----------------------------+
| 13-16            | Echogram         | Presently unused           |
+------------------+------------------+----------------------------+

Care should be taken to preserve the FID in all of the files that follow,
because they will be important for matching the pieces of data properly.

UTM coordinates
~~~~~~~~~~~~~~~

As can be seen by viewing the data in one of the ways above, the radar data
contains geographical longitude and latitude coordinates. It's frequently
easier to work on a projected coordinate system, so I would run

::

    h5_add_utm.py data/gl3_radar_2012.h5 data/gl3_radar_2012_utm.h5

to create a copy of the file with UTM coordinates appended. Repeating the
previous step, it can be seen that the *eastings*, *northings*, and *zones*
columns are now populated. **This assumes that you've installed pyproj**
(see Dependencies_).

-  Note: ``h5_add_utm`` is pretty dumb about UTM zones. It works fine most
   places, but if you're working in SW Norway or Svalbard, there are exceptions
   to the normal 6\ :math:`^\circ` grid and you might need to tweak the code.

Pre-caching (optional)
~~~~~~~~~~~~~~~~~~~~~~

Finally, as an optional step, I'm going to generate a cached copy of the radar
data. This will speed things up while picking, and allows me to do some initial
preprocessing to remove bad radar soundings, etc.

::

    h5_generate_caches -d cache -g -b -r --dc=0 data/gl3_radar_utm.h5
    h5_generate_caches -d cache -g -b -r --dc=1 data/gl3_radar_utm.h5

This creates a preprocessed binary copy of each line in the directory ``cache``
for which

-  fixes for when the GPS was in "car-mode" (static GPS) have been attempted
-  repeated soundings in the same location have been filtered out
-  blank traces from caused by triggering errors have been removed

The ``--dc`` switch specifies which channel to operate on, which is only
important for the "dualdar" system (otherwise it should be ``0``, which is the
default).

To see what preprocessing options are available, type ``h5_generate_caches``
without any argument.

Ice thickness picking
---------------------

Basics
~~~~~~

Picking is the process of measuring the time before the arrival of each
interesting return wavelet. Picking is labour-intensive, although I tried to
automate the easy parts.

To get started, run

::

    icepick -f data/gl3_radar_utm.h5

The way this works is that the window that opens shows a grey-scale *radargram*
in the top panel, with the eight individual traces in the bottom panel. The
location that the traces are from is shown by the vertical yellow lines in the
radargram. Assuming there are more than eight traces (normally the case), the
display can be panned across the radargram with the **h** and **l** (ell) keys.
In my case with the 2012 radar data from Glacier 3, the first line only
contains a single trace, so panning doesn't do anything.

The terminal in which ``icepick`` was launched now accepts icepick-specific
commands. Typing

::

    info

gives information about the current line. For this data, it tells me

::

    data/gl3_radar_2012_utm.h5
    line: 0
    # traces: 1
    # samples: 256
    sample interval: 4e-09 s
    depth resolution: 0.336 m
    vertical range: 86.016 m
    pick-mode: bed

From top to bottom, this tells me what file I'm operating on, the line number
(starts at 0, as in the HDF dataset), the number of traces (``nx``), the number
of samples per trace (``nz``), the sampling interval, and estimates of the
vertical resolution and the maximum depth imaged, assuming the material is ice.
The final line, ``pick-mode``, indicates that any picks we perform now are for
the glacier ``bed`` (more on that in a moment).

Typing

::

    help

gives a (potentially non-exhaustive) list of valid commands. To switch to line
#1, type

::

    open 1

The process of picking is fairly simple. In the lower panel of the icepick
window (where the individual traces are shown), click the mouse on the part of
the trace representing a reflection to be timed. Right-clicking removes the
pick if you've made a mistake. Fine adjustments can be made by pressing the
**j** (down) and **k** (up) keys. Whenever the side-scrolling keys are pressed
(**h** and **l**), a line representing the picks is drawn on the radargram.
Presumably, the bed should be picked on every trace where it can be identified.

Once the picks are satisfactory, type ``save`` to save the timing data to the
folder ``picking``. If it spits back something like

::

    IOError: [Errno 2] No such file or directory:
    'picking/gl3_radar_2012_utm_line1.csv',

then you probably forgot to create the ``picking`` folder in `the section on
set-up <#set-up>`__.

Filtering
~~~~~~~~~

There are a number of filters that can be applied with the ``f`` command, using
the syntax

::

    f FILTERNAME

Some common filter names are:

-  ``dewow``: applies a "dewowing" highpass filter
-  ``lowpass``: applies a generic frequency lowpass filter
-  ``lowpass_ma``: applies a generic time-domain lowpass filter
-  ``gc``: applies a linear gain control
-  ``agc``: applied a nonlinear automatic gain control (usually more fun than
   useful)
-  ``fkmig``: performs F-K (Stolt) migration, and takes a sample number as an
   optional argument indicating time zero (the airwave)

Furthermore,

-  Typing ``f`` without an option lists the filter history, so you can see
   exactly how the current data has been modified.

-  Typing ``nf`` undoes all filter effects (except for those that happened
   during cache-generation or automatically when loading the line), and
   restores the original data.

There are lots of other filters. All filters are defined in the file
``filter_defs.py``, which is in the place where ``irlib`` is installed.
Modifying this file permits custom filters to be defined.

A final adjustment is ``gain``, which adjusts the display contrast of the
radargram. All filters accessed through ``f`` or ``gain`` are reversible, so
there is no risk of permanently damaging the data by experimenting.

For bed picking, **I strongly recommend performing some kind of waveform
migration (``fkmig`` is a good option).** Migration attempts to restore
reflector geometries, and is important wherever bed slopes may be large. The
gain control filters combat signal attenuation, and are also frequently useful.
The various bandpass filters are selectively worthwhile, but may introduce
artefacts in the reflection and direct-coupling wavelets, so some caution is
required.

Direct coupling
~~~~~~~~~~~~~~~

In order for timing data to be generated, a reference time must be known.
Because it's not easy for us to know the exact time that the transmitter
emitted a pulse into the ice, we use the airwave as a timing reference. The
airwave travels directly from the transmitting antennas to the receiving
antennas at the speed of light (:math:`\approx 3\times10^8\text{ m\,s}^{-1}`,
so the emission time can be calculated by knowing the airwave arrival time.

To switch to direct-coupling mode, type

::

    mode dc

and a label should appear in picking window indicating the mode change. All
picks made in ``dc`` mode will have a red dot rather than blue.

To change back to bed mode, type ``mode bed``.

Automated picking
~~~~~~~~~~~~~~~~~

To save time, picking can be done automatically. For example, to automatically
pick the airwave across the whole radar line, use the ``autodc`` command. If me
know that the airwave is between samples 75 and 125 (right vertical axis on the
radargram), then we can give this as a hint by typing

::

    autodc 75 125

*icepick* then uses a set of heuristics to try and figure out where the airwave
is in each trace, subject to the vertical constraints.

-  There is a minimum vertical range for the algorithm to work. I forget what
   it is, but it's something around 20. If ``autodc`` doesn't work, try
   increasing the range arguments.

Automatically picking the airwave usually works pretty well. Automatically
picking the bed reflection is more hit-and-miss. The command ``autobed`` works
pretty much the same way as above, and usually does a decent job when the
radargram is very clear. Even when the radargram is more complicated, I usually
give ``autobed`` a shot, and then go through making the (many) necessary
corrections.

Pick rating
-----------

Rating is used to quantify the certainty of each pick. I use the following
rating table

+----------+---------------------+
| Rating   | Approximate Error   |
+==========+=====================+
| 5        | 1.4 m               |
+----------+---------------------+
| 4        | 1.7 m               |
+----------+---------------------+
| 3        | 2.2 m               |
+----------+---------------------+
| 2        | 3.5 m               |
+----------+---------------------+
| 1        | 7.1 m               |
+----------+---------------------+

Ratings could be tabulated manually. For efficiency, I use a program similar to
*icepick*

::

    icerate -f data/gl3_radar_2012_utm.h5

but this program is not polished to the same standard as *icepick* and
*irview*.

Ice thickness calculation
-------------------------

Antenna spacing
~~~~~~~~~~~~~~~

A last ingredient before ice thickness can be calculated is an *offsets* file,
which contains information about how much antenna spacing there was for each
line. Hopefully this information is contained in field notes. Open the CSV
created previously with ``h5_dumpmeta``. Save it as XXXX\_offsets.csv, where
XXXX is the prefix of the HDF file (also prepended onto the picking and rating
files created above). In this case, it should be saved as
``gl3_radar_2012_utm_offsets.csv``. Keep the first column ("FID") and delete
all other columns. Add a new column, and add the appropriate antenna spacing in
meters for each row (see `FID interpretation <#fid-interpretation>`__). Delete
the header row, save, and exit. Convert the CSV to a Tab-delimited file, e.g.:

::

    cat gl3_radar_2012_utm_offsets.csv | sed 's/,/\t/g' > \
        gl3_radar_2012_utm_offsets.txt

Data join
~~~~~~~~~

Calculating ice thickness is fairly trivial, so the only challenge is in
properly integrating all of the data. The steps are:

-  Take all soundings for which both a pick and a rating exist
-  Find the proper antenna spacing
-  Assuming an ice velocity, calculate reflector depth with the Pythagorean
   theorem

I use the script in ``scripts/radar/join_radar.py`` to do all of this.

::

    python join_radar.py gl3_radar_2012_utm data/gl3_radar_2012_utm.h5

which should generate a file containing data similar to:

+-----------+------------+-------------+------------------+
| easting   | northing   | depth (m)   | variance (m^2)   |
+===========+============+=============+==================+
| 609481    | 6760344    | 25.29       | 3.125            |
+-----------+------------+-------------+------------------+
| 609477    | 6760340    | 23.18       | 3.125            |
+-----------+------------+-------------+------------------+
| 609473    | 6760339    | 23.88       | 3.125            |
+-----------+------------+-------------+------------------+
| 609470    | 6760337    | 24.59       | 3.125            |
+-----------+------------+-------------+------------------+
| 609467    | 6760336    | 24.59       | 3.125            |
+-----------+------------+-------------+------------------+
| 609463    | 6760335    | 26.69       | 5.55555555556    |
+-----------+------------+-------------+------------------+
| 609458    | 6760333    | 26.69       | 5.55555555556    |
+-----------+------------+-------------+------------------+
| 609449    | 6760329    | 29.47       | 5.55555555556    |
+-----------+------------+-------------+------------------+
| 609444    | 6760327    | 30.85       | 5.55555555556    |
+-----------+------------+-------------+------------------+
| 609441    | 6760324    | 34.3        | 3.125            |
+-----------+------------+-------------+------------------+
| 609437    | 6760323    | 34.99       | 3.125            |
+-----------+------------+-------------+------------------+
| 609434    | 6760322    | 32.92       | 3.125            |
+-----------+------------+-------------+------------------+

Raster interpolation
--------------------

The general interpolation scheme is discussed in `Interpolation
<#interpolation>`__. A brief description and the commands I used to generate a
bed map are given below.

Mask file
~~~~~~~~~

I generate a mask covering the area of Glacier 3 based on the outline traced
from satellite imagery. This provides a domain for the interpolation scheme.
Using the outline shapefile from `Outlines <#outlines>`__:

::

    gdal_rasterize -of GTiff -a id -tr 20 20 -te 606100 6757400 611100 6760800\
                   -l outline_gl3 outline_gl3.shp mask_gl3.tif
    gdal_translate -of AAIGrid mask_gl3.tif mask_gl3.asc

Data concatenation
~~~~~~~~~~~~~~~~~~

Since ice thickness needs to be zero at the glacier margin (assuming no cliffs
or steep bulges), I append the depth sounding data generated `above
<#data-join>`__ with samples taken from the glacier margin. I produced the
margin file using a GIS, and prescribed a depth of 0 m and a variance of 0.1 m
at every point (*gstat* doesn't like zero uncertainties). Then,

::

    cat depth_gl3_radar_2012_utm.xyz gl3_outline_100m.xy > \
        kriging/gl3_depth_outline_2012.xyz

Variogram estimation
~~~~~~~~~~~~~~~~~~~~

I created a proto-\ *gstat* configuration file called ``gl3_12_2p.gst`` and
containing the lines:

::

    data(gl3): 'depth_outline_2012.xyz', x=1, y=2, v=3, V=4, d=2, \
               average=1, max=100, radius=1000;
    set zero=20;

The first line creates a datasource from the concatenated ice thicknesses, and
indicates that the columns correspond to *x* and *y* spatial coordinates, the
interpolated value (*v*), and the variance (*V*), respectively. The argument
``d=2`` assumes a quadratic trend, ``average=1`` permits averaging of points
that are very close, ``max=100`` sets a maximum number of observations for each
interpolated point, and ``radius=1000`` sets a maximum search neighbourhood.

The second line declares that points within 20 m are indistinguishable from
each other.

Running this

::

    gstat gl3_12_2p.gst

opens an interactive *gstat* session, from which variogram estimates can be
saved. I assume that, because Glacier 3 is roughly east-west oriented, the
variogram should be split into east-west and north south components, and I save
a variogram estimate for each.

Model variogram fitting
~~~~~~~~~~~~~~~~~~~~~~~

Variogram fitting can be performed in *gstat*, but I use a Python script
(``fit_variogram.py``) because it gives me more control over the fitting
routine and is more suited for anisotropic variograms than the built-in tools.

.. figure:: images/variograms.png
   :alt: Estimated experimental variograms (points) and modelled
   variograms (lines) for the major and minor axes of Glacier 3

   Estimated experimental variograms (points) and modelled variograms
   (lines) for the major and minor axes of Glacier 3

Once a suitable model variogram has been found, the *gstat* configuration file
can be modified:

::

    variogram(gl3): 1439 Sph(1043.8, 90, 0.4428) + 514 Sph(151.7, 0, 0.9750);
    mask: 'mask_gl3.asc';
    predictions(gl3): 'predictions/pred_gl3_12.asc';
    variances(gl3): 'variances/var_gl3_12.asc';

Running this again will perform the interpolation. See the *gstat* manual for
details.
