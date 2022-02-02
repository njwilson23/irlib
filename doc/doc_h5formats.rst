Changes in IEI h5 file format
=============================

The following is a brief overview of the various different hdf file formats that have come about from
different versions of the IEI Software for acquiring BlueSystem Radar data. 

These modifications have broken irlib code in the past but version 0.5 works-around these idiosyncrasies. 

An easy way to examine the format of the h5 file you have is to look at it with an hdf viewer like vitables. 


Version 5.1 IEI (Sept 2016)

	Added correction for lat and lon such that the western and southern hemispheres are negative numbers
	Lat and Long are the field names and the values are floating point numbers

	PCSavetimestamp attributes contains the same string as before but the timestamp is before the 
	GPSCaptureEvent_StartBufferCapture.ms field 	The format is dd/mm/yyy_hh:mm:ss

Version 4 IEI: Format of PCSavetimestamp string is:

	With no GPS used: GPSCaptureEvent_StartBufferCapture.ms:-99,BufferCaptureTime.ms:264,PPS_NO
	With standard GPS reading: GPSCaptureEvent_StartBufferCapture.ms:72,BufferCaptureTime.ms:336,PPS_NO
	With PPS GPS reading: GPSCaptureEvent_StartBufferCapture.ms:72,BufferCaptureTime.ms:336,PPS_YES
    
	Importantly, in some cases the actual timestamp is saved as a comment field. This is difficult to get!
	The format is dd/mm/yyy 

Version ? IEI: (ca. 2012?)

	lat and lon were both recorded as positive numbers in fields Lat_N and Long_W. If either the lat and 
	lon are negative numbers, then irlib assumes that this paradigm is not in effect. 
	The user can use --swap_lon and --swap_lat to change the sign of either lat and lon in h5_add_utm.py 
	and --positivecoords in h5_replace_gps.py
	
	Computer time stored as PCSavetimestamp: mm/dd/yyyy_hh:mm:ss like "3/12/2014_11:49:20 AM" There are no
	GPSCapture stats.
	
