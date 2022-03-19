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

::

    SYNTAX: h5_consolidate INFILE1 INFILE2 [...] -o OUTFILE

		Combines multiple datasets (>1) into a single concatenated dataset.

``h5_consolidate`` combines multiple datasets into a single dataset. In the
process, lines are re-numbered so that they stay in sequential order.
Concatenating datasets is useful, for example, to combine multiple surveys
collected on different days into a single file that is easier to manage (but
larger).

h5_replace_gps
~~~~~~~~~~~~~~

::

    SYNTAX: h5_replace_gps infile outfile gpsfile {gpx,ppp} {iprgps,iprpc,both} [OPTIONS]
	
		This tool replaces the existing geographical data in a ice radar HDF
		database with data taken from a GPX file, e.g. obtained from a handheld or
		external GPS unit or from a CSV file, e.g. obtained from a PPP output of 
		GPS data

	Positional arguments:
		infile           	input HDF (.h5) filename, with or without path, 
							for which GPS or PC timestamps exist
		outfile          	output HDF (.h5) filename, with or without path, 
		                    if this file exists, it will be overwritten
		gpsfile          	GPS filename(s), with enhanced location, with or 
		                    without path / wildcards
		{gpx,ppp}        	Select which format the gps file is in - either 
		                    gpx or ppp
		{iprgps,iprpc,both}	Select which timestamp to match gps timestamps to 
		                    - iprgps (recommended), iprpc (if iprgps not available) 
							or both (use caution)

    Optional arguments:
		-t hh 	The hour offset (hh) of the GPR computer from UTC (default = 0)
		-l n    Work only on line (n); default works on all lines
		-d n 	Set the max time delta permissible for matching locations to (n) 
				seconds; default is 15 seconds
		-o n 	Adds an offset (n) to the elevations to account for the height of 
				GPS off the ice or different geoid, use a neg. number to subtract.
		-n  	Replace coordinates in HDF with no appropriate supplementary GPS 
				counterpart with 'NaN'. By default, the original coordinates are retained.
		-p  	Keep all coordinates positive (use with old h5 format where 
				Lat_N and Long_W).
		
If GPS data collected from the on-board receiver are missing or of poor
quality, they can be replaced by data from a hand-held GPS receiver. The data
from the hand-held receiver must be exported as or converted to GPX format,
which is a standard open format. Calling ``h5_replace_gps`` creates a copy of
the original dataset with the new coordinates inserted. Command-line flags can
be used to specify matching tolerances and which lines to work on.

h5_add_utm
~~~~~~~~~~

::

    SYNTAX: h5_add_utm INFILE OUTFILE

        Replaces geographical coordinates in INFILE with UTM coordinates
        in OUTFILE. Does not perform any datum shift. Projection is calculated
        assuming that the data from neither from western Norway nor Svalbard.

``h5_add_utm`` uses the *pyproj* library to append projected UTM zone
coordinates to datasets that only include lon-lat coordinates. This is a
required step for many of the data processing operations that might be used
later.

The UTM zone is calculated based on a naive algorithm that is ignorant of the
exceptional UTM circumstances in the vicinity of western Norway and Svalbard.

Works with 2 formats from BSI HDF files: 
  	Old format - Latitude and longitude data in BSI HDF files are unsigned. It 
		is assumed to be in the western hemisphere by default. Passing the --swap_lon 
		key forces longitudes to be interpretted from the eastern hemisphere.
		UTM projection is calculated assuming that the data from neither from western 
		Norway nor Svalbard.
	New format - Latitude and longitude data in BSI HDF files are signed to indicate 
		hemisphere. If any lat or lon values are negative, the --swap_lon key is disabled

h5_generate_caches
~~~~~~~~~~~~~~~~~~

::

    SYNTAX: h5_generate_caches HDF_SURVEY [OPTIONS]

        -d [DIR]    cache directory (default: cache/)
        -g          fix static GPS issues
        -s          smoothen coordinates
        -b          remove blank traces caused by triggering failure
        -r          remove stationary traces by averaging all traces within # m 
					(defaults to 0 m or off), recommend 3 for L1 GPS
        -f          force regeneration of existing caches
        -q          silence standard output
        -e          print failed datacaptures
        --dc=[#]    specify datacapture (default: 0)
		-n 			remove traces with NaN coordinates
		-i			interpolate over NaN coordinates (overrides -n)
		-v			print failed datacaptures

Caching improves performance and is a very good idea. ``h5_generate_caches``
creates caches (``.ird`` files) for every line within a survey, and optionally
applies a number of pre-processing steps to the data:

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

h5_dumpmeta
~~~~~~~~~~~

::

    SYNTAX: h5_dumpmeta infile [OPTIONS]

    Positional arguments:
		infile	input HDF (*.h5) filename, with or without path, if you use 
		wildcards 
				in linux, put this in quotes

    Optional arguments:
		-o 		output file BASENAME [if missing, will be automatically 
				generated]
		-c 		create csv metadata file
		-w 		create a waypoint metadata shapefile
		-l 		create a line metadata shapefile
		--clobber  	overwrite existing files
		

``h5_dumpmeta`` exports the radar metadata to a CSV file. The actual sounding
data is not included.

h52mat
~~~~~~

::

    SYNTAX: h52mat SURVEYFILE OUTFILE [options]

    SURVEYFILE is the HDF5 file generated by IceRadar.
    OUTFILE is the anme of the *.mat file to be generated.

    Options:
        g       fix static GPS issues
        s       smoothen coordinates
        b       remove blank traces (trigger failure)
        r       remove stationary traces
        o       overwrite
        q       silence standard output

``h52mat`` converts HDF data to a MATLAB ``.mat`` file. The filters from
``h5_generate_caches`` are available. For those who prefer MATLAB, the rest of
this document can be ignored.


Recommended data cleaning workflow
----------------------------------

The following steps are very helpful for data cleaning and streamlining
workflow. Also some of the steps are prerequisites for subsequent
analyses, so **do this in the correct order**. It is really very
important that you take notes on what you did so that your workflow
can be recreated later. It is recommmended you open a document and copy paste
what you did from the terminal in there for safekeeping. Also, you can
copy the screen output there too. As you go be aware that some scripts
will overwrite files. Recommend that you use unique file names that
represent the step that you just completed.

-  ``h5_dumpmeta.py``
-  ``h5_consolidate``
-  ``h5_replace_gps.py``
-  ``h5_add_utm.py``
-  ``h5_dumpmeta.py``
-  ``h5_dumpmeta.py``

Once this has been been completed the data is ready to be used in IcePick2.py, 
which will be elaborated on in the next chapter.

