Command-line Utilities
======================

The command-line utilities in *radar_tools* are useful for performing data
management and pre-processing tasks on HDF radar datasets, as well as for
performing basic data exploration and conversion tasks.

In general, typing any of the utilities without arguments yield invocation and
usage instructions that are printed to the screen. This section summarizes
individual tool's functionality.

Data management
----------------

h5_consolidate
~~~~~~~~~~~~~~

``h5_consolidate`` combines multiple datasets into a single dataset. In the
process, lines are re-numbered so that they stay in sequential order.
Concatenating datasets is useful, for example, to combine multiple surveys
collected on different days into a single file that is easier to manage (but
larger).

h5_replace_gps
~~~~~~~~~~~~~~

If GPS data collected from the on-board receiver are missing or of poor
quality, they can be replaced by data from a hand-held GPS receiver. The data
from the hand-held receiver must be exported as or converted to GPX format,
which is a standard open format. Calling ``h5_replace_gps`` creates a copy of
the original dataset with the new coordinates inserted. Command-line flags can
be used to specify matching tolerances and which lines to work on.

h5_add_utm
~~~~~~~~~~

``h5_add_utm`` uses the *pyproj* library to append projected UTM zone
coordinates to datasets that only include lon-lat coordinates. This is a
required step for many of the data processing operations that might be used
later.

The UTM zone is calculated based on a naive algorithm that is ignorant of the
exceptional UTM circumstances in the vicinity of western Norway and Svalbard.

h5_generate_caches
~~~~~~~~~~~~~~~~~~

Caching (discussed later?) improves performance and is a very good idea.
``h5_generate_caches`` creates caches (``.ird`` files) for every line within a survey, and
optionally applies a number of pre-processing steps to the data:

    - **static gps correction**: attempt to recognize period when the GPS was
      in "static mode", and interpolate continuous positions.

    - **smoothen coordinates**: filter noisy position data

    - **remove blank traces**: exclude empty soundings from the cache

    - **remove stationary traces**: attempt to recognize period when the radar
      sled was motionless, and remove redundant soundings

``h5_generate_caches`` should be the last of the data management scripts to
run, because modifying the original HDF dataset won't affect the caches until
they are regenerated.


Exploration and conversion
---------------------------

h5_dump_meta
~~~~~~~~~~~~

``h5_dump_meta`` exports the radar metadata to a CSV file. The actual sounding
data is not included.

h52mat
~~~~~~

``h52mat`` converts HDF data to a MATLAB ``.mat`` file. The filters from
``h5_generate_caches`` are available. For those who prefer MATLAB, the rest of
this document can be ignored.

