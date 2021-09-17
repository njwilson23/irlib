.. _data_cleaning_and_workflow:

Data cleaning and workflow
--------------------------

The following steps are very helpful for data cleaning and streamlining
workflow. Also some of the steps are prerequisites for subsequent
analyses, so **do this in the correct order**. Note that these commands
are intended for the most recent irlib by Derek. It is really very
important that you **take notes on what you did** so that your workflow
can be recreated later. Recommmend you open a document and copy paste
what you did from the terminal in there for safekeeping. Also, you can
copy the screen output there too. As you go be aware that some scripts
will overwrite files. Recommend that you use unique file names that
represent the step that you just completed.

The picking files (see picking subfolder) represent your hard work.
**Make sure you save your pick files carefully (versioned zipfile?).**

-  *Copy* original h5 files into a data directory where you have space,
   write permissions and can back up.
-  *Copy* ppp csv files to a subdirectory called ppp
-  Create a folder called *metadata* (recommended to keep things tidy)
-  Several subfolders will be created by irlib eventually:

   -  cache
   -  picking
   -  results

-  Depending on what steps you take you may also see the following
   subdirectories

   -  rating
   -  offsets

To view and manually edit HDF5 files use **ViTables**, simply open the
GUI using the command below. Once open, choose the file you would like
to browse. This is not necessary for data cleaning, but is helpful for
viewing HDF5 files.

``vitables``

**All the following programs have settings that you can find out about
in the irlib manual or by typing -h after the program in the terminal.
Only the most recommended settings are provided below**

**h5_dumpmeta.py** creates a csv of all the metadata available. It will
output a sample of the data to the screen type -c to save a csv version,
-w to save a point shapefile and -l to save a line shapefile

``h5_dumpmeta.py "raw_ipr*.h5" -c -w``

**h5_consolidate.py** takes a bunch of h5 files and concatenates them
together. This means new line numbers will be assigned (**It is highly
recommended that you note down the files that you consolidate and their
respective line numbers**). Also a very large h5 file will be generated.
(note wildcards for filenames is ok)

``h5_consolidate.py "raw_ipr*.h5" -o ipr_survey.h5``

- **h5_replace_gps.py** takes a gpx trk file or a `NRCan PPP`_ file and
replaces the onboard gps data. Run this if the onboard GPS was not on or
you have ancillary position data that is better quality than the
original onboard gps data. This creates new file.

-  **-o** accounts for offset from the gps to the ice (depends on ipr
   set up), this value should be negative (example here the GPS was 18
   cm above the ice surface so subtract 0.18 to get the elevation of the
   ice).
-  **-n** replaces the location data that doesn't match by timestamp to
   the supplementary gps data with NAN -- to be conservative choose this
   option
-  **-d** the delta time that constitutes a match, can be changed from
   the 15 second default (If you are using 1 Hz Topcon data you should
   reduce this to 1-2 seconds at the most). Check the output and how
   many traces were modified to see if that worked ok.
-  **timesource** select which time stamp to match gps timestamps to:
   iprgps (recommended - will match time to the internal IPR GPS
   timestamp), iprpc (matches to the IPR PC clock -- only if iprgps is
   NOT available), or both (use caution since this may cause the track
   to go back and forth -- not recommended).

The -d and -o values in this command are a common example, however,
there are many forks in the road.

``h5_replace_gps.py ipr_survey.h5 ipr_survey_ppp.h5 "nrcan*.csv" ppp iprgps -d 1 -n -o -0.18 ``

**h5_add_utm.py** converts lat/lon to UTM Cartesian coords so they are
easier to deal with. This creates new file. You must do this to go on
further.

``h5_add_utm.py ipr_survey_ppp.h5 ipr_survey_ppp_utm.h5``

**h5_dumpmeta.py** Dump metadata again to look at the timestamps or
check field notes. Confirm all is well by looking at the output of this
script (compare to the earlier version of dumpmeta).

``h5_dumpmeta.py ipr_survey_ppp_utm.h5 -o metadata/ipr_survey_ppp_utm -w -c``

**h5_generate_cache.py** creates caches (python pickle \*.ird files) of
each line and performs some further cleaning steps.

-  **-r** removes traces that are within a given distance by averaging
   them together. If you have low precision GPS that number can be 3-5 m
   but if you want to get high resolution data and you are using a good
   GPS, that number can be lower.
-  **-n** removes traces with NaN coordinates - to be conservative,
   choose this.
-  **-i** interpolates NaN coordinates
-  **-s** smooth coordinates - I don't recommend this
-  **-f** keep this on if you have a cache already.
-  other switches are not required

``h5_generate_caches.py ipr_survey_ppp_utm.h5 -r 1.0 -n -b -f``

**icepick2.py** This script handles the picking and display of data.

-  Open each line in icepick and pick the DC and Bed for all traces that
   you are confident in
-  See below for how to do this
-  Make sure to save the picking subfolder so you can preserve your work


**join_radar.py** The final step is to merge the picking data with the
radar metadata to give the ice thickness. If the entire survey (all
lines) have the same radar offset (separation distance) then you can add
this here. In the example below, we'll use the default radar velocity
and a separation of 15 m. Also the -n will remove any traces that have
no pick data. If you want to see areas that were not picked don't use
this. The -o value in this command is a common example, however, there
are many forks in the road.

``join_radar.py ipr_survey_ppp_utm.h5 -c -w -o 15 -n``

You will see the output in the results folder.

.. _picking_data:

Picking Data
------------

Picking is time-consuming and challenging. Use icepick for an
interactive experience with the data.

``icepick2.py -f mydata.h5``

-  Type help for help on available functions
-  gain - this is a contrast adjustment on the display
-  filter **help filtertype** to get more info...
-  type **filter** (or f) without listing a filtertype to see what's
   applied
-  type **nf** to remove all filters
-  remember filters act one on top of the other, so order will be
   important.
-  you can click on the radargram to see the fid, sample and time (fid
   is Nat's unique ID for each trace)

-  type map on (off) to open (close) map window
-  type pick on (off) to open (close) pick window
-  pick bed with mouse

   -  middle mouse button switches mode to-from DC/bed
   -  move pick with j and k one at a time / fine-tune
   -  h and l move next set of traces

-  pick save – adds folder picking and a file with fid, start, end and
   difference (twtt in # of samples)
-  pick load – will load a pick file

-  Auto pick
-  **pick dc** will pick top automatically ('first break' method)
-  **pick bed** use mouse to figure max and min sample number off the
   radargram
-  **pick bed 33 89** - picks bed between these two sample numbers
-  fix picks after by deleting with with the right mouse button

.. _suggested_picking_workflow:

Suggested picking workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Open each line in Icepick2
-  Pick DC - if this can't be done auto then there is a problem.
   (ideally do this pick *before* filtering b/c filters change data)
-  Fool around with filters (find the best combination for optimal
   visibility - record this setting for each line - or portion thereof)

   -  try dewow
   -  try gc 1 to 2 - usually fine
   -  try ringing filter (but not sure about this)
   -  try slight lowpass_td (like 7) since that is good sometimes before
      autopicking (but this may add a tiny wave to DC - see above)

-  Next pick the bed

   -  Try to automate this and then fix the picks as you go


