Changes in IEI h5 file format
=============================

The following is an overview of the various different hdf file formats that have come about from
different versions of the IEI Software

timestamps taken as ddmm (newer) or mmdd format (older) and metadata outputs in ISO format


 (as per Laurent Mingo)

Version 5.1 IEI (Sept 2016)
- Added correction for Lat Lon,
    E and W with positive and negative signs

Version 4.4.1 IEI: Format of GPS string is:

    With no GPS used: GPSCaptureEvent_StartBufferCapture.ms:-99,BufferCaptureTime.ms:264,PPS_NO
    
    With standard GPS reading: GPSCaptureEvent_StartBufferCapture.ms:72,BufferCaptureTime.ms:336,PPS_NO
    
    With PPS GPS reading: GPSCaptureEvent_StartBufferCapture.ms:72,BufferCaptureTime.ms:336,PPS_YES

lat lon format -- (which as of IEI Version 5.1, Sept 2016 is a signed float)
