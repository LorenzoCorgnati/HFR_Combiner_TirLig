# HFR_Combiner_TirLig

This application is written in Matlab language and it is based on HFR_Progs_2_1_2 and M_Map toolboxes.
This application is designed for High Frequency Radar (HFR) data management according to CNR-ISMAR processing workflow.

HFR_Combiner_TirLig_v30.m is the main code file.

This application automatically reads Radial data from HF radar
installation in Tyrrhenian and Ligurian Sea (MONT, TINO), generates totals
and stores radial and total files.
Grid: 2 km - Combination Radius: 3 km.

The 1.1 release is able to recover old radial data used for uncomplete
total maps generation in order to remake totals using new radial data
sent with a delay by the communication system.

The 1.1.x release works on shared folder radardisk and radarserver and
plots total maps in an optimized way (reduced area and more visible
arrows) and saves incremental copies of RadStruct.mat file.

The 1.1.5 release reads the input parameters from an external text file,
namely the folder paths, the latitude and longitudes values for maps
representation and for total maps evaluation, the velocity threshold for
totals cleaning, the spatial threshold for totals generataion and the
time period limits for totals evaluation.

The 1.2 release starts from the architecture of 1.1.5 version and
implements cleaning and corrections on the radial data in order to
operate a first level Quality control.

The 1.3 release stores hourly total data in netCDF format.

The 1.3.1 release stores hourly total data in netCDF format on
the THREDDS server, into the RITMARE dedicated folders.

The 1.3.2 release runs a function at its startup in order to delete from
the RadarDisk RadialRealTime folder the radial files older than a week
(Archivalist rsyncs radial files of the last week). This task speeds up
the directory lisiting and RadStruct operations.

The 1.3.4 release generates total velocity files with current
measurement unit in m/s according to the EmodNet requirements.

This release stores the tuv and the netCDF total files in a folder structure
with the scheme yyyy/yyyy_mm/yyyy_mm_dd/ locally and on the RadarDisk,
and with the structure Last/yyyy_mm_dd/ on the THREDDS server.

The 1.3.5 release makes local copies of the radial files locally and
processes them locally. Then the resulting total products are moved in
the destination folders. The local copies are deleted after processing.

The 2.0 release creates netCDF 4 files containing total velocities, U and V
errors, UV covariance and GDOP. The file format is compliant with CF-1.6
Unidata Data Discoverty, INSPIRE and ACDD conventions.

The 2.1 release stores hourly radial and total data in netCDF format on
the THREDDS server, into the RITMARE dedicated folders.

The 2.2 release stores hourly radial and total data in netCDF format in a
folder structure compliant to the ISMAR THREDDS requirements.

The 3.0 release implements the data and metadata structure compliant to
CMEMS needs and performs QC tests.

Author: Lorenzo Corgnati
Date: January 26, 2017

E-mail: lorenzo.corgnati@sp.ismar.cnr.it
