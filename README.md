## Sun Nuclear Data Extraction Tools 

by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2015, University of Wisconsin Board of Regents

The Sun Nuclear Data Extraction Tools are a compilation of functions that extract data from [Sun Nuclear](http://www.sunnuclear.com) application and detection system exported files into MATLAB structures and arrays, including the SNC ArcCHECK&reg; and IC Profiler&#8482; systems.  These tools are used in various applications, including [viewray_mlc](https://github.com/mwgeurts/viewray_mlc), [viewray_fielduniformity](https://github.com/mwgeurts/viewray_fielduniformity), and [viewray_radiso](https://github.com/mwgeurts/viewray_radiso).

ArcCHECK and IC Profiler are trademarks of Sun Nuclear Corporation.

## Contents

* [Installation and Use](README.md#installation-and-use)
* [Compatibility and Requirements](README.md#compatibility-and-requirements)
* [Tools and Examples](README.md#tools-and-examples)
  * [ParseSNCacm](README.md#ParseSNCacm)
  * [ParseSNCprm](README.md#ParseSNCprm)
  * [ParseSNCtxt](README.md#ParseSNCtxt)
  * [AnalyzeProfilerFields](README.md#AnalyzeProfilerFields)
  * [AnalyzeACFields](README.md#AnalyzeACFields)
* [Event Calling](README.md#event-calling)
* [Gamma Computation Methods](README.md#gamma-computation-methods)
 
## Installation and Use

To install the Sun Nuclear Data Extraction Tools, copy all MATLAB .m files from this repository into your MATLAB path.  If installing as a submodule into another git repository, execute `git submodule add https://github.com/mwgeurts/snc_extract`.

## Compatibility and Requirements

The Sun Nuclear Data Extraction Tools have been validated with SNC Patient version 6.2.3 and SNC Profiler version 3.3.1 using MATLAB versions 8.3 through 8.5 on Macintosh OSX 10.8 (Mountain Lion) through 10.10 (Yosemite).  These tools do not require any non-standard MATLAB/Java libraries or functions.

All functions are compatible with both JVM and non-JVM (MATLAB was executed with the `-nodisplay`, `-nodesktop`, or `-noFigureWindows` flags) environments.  If a display is present, certain functions will display a progress bar while executing.

## Tools and Examples

The following subsections describe what inputs and return variables are used, and provides examples for basic operation of each tool. For more information, refer to the documentation within the source code.

### ParseSNCacm

`ParseSNCacm()` extracts data from a SNC ArcCHECK movie file (.acm) and returns the contents of the file as a MATLAB structure.  See below for a full list of the structure fields returned.  This function will display aprogress bar while it loads (unless MATLAB was executed with the `-nodisplay`, `-nodesktop`, or `-noFigureWindows` flags).

The second input parameter can either contain one or multiple files.  If multiple files are selected, the returned data is a concatenation of all measured profiles.  However, only the header data (array calibration, background, ignore flags, etc) for the last file is returned.

If a given field is not found, this function will gracefully ingore it and the returned structure will not contain the field. If the field is found but no contents were specified, the returned field will be an empty cell array.

This function has been tested with SNC Patient version 6.2.3 and .acm file rev C. Note, there are additional fields in SNC .acm files that are not currently imported by this function. Additional fields can be added using the modular `search` variable, declared within this function. Refer to the documentation within the source code for more information.

The following variables are required for proper execution:

*  path: string containing the path to the DICOM files
* names: string or cell array of strings containing the file(s) to be loaded

The following structure fields are returned upon successful completion:

* filerev: string containing the file revision letter
* filename: string containing the original filename
* timestamp: date and time .acm file was saved, as an integer
* institution: string containing the institution
* dcal: string containing the calibration file name
* version: string containing the SNC software version
* dorientation: string containing the ArcCHECK orientation
* dssd: SSD, in cm
* shiftx: phantom X shift, in mm
* shifty: phantom Y shift, in mm
* shiftz: phantom Z shift, in mm
* rotx: phantom X rotation, in deg
* roty: phantom Y rotation, in deg
* rotz: phantom Z rotation, in deg
* dmodel: string containing the detector model
* dserial: string containing the detector serial number
* dfirmware: string containing the detector firmware revision number
* dinterval: collection interval, in ms
* dthreshold: background threshold, in percent
* dplug: string indicating whether plug was present for acquisition
* dplugdose: measured dose (in Gy) to central cavity, if provided
* mroom: string containing the treatment room
* mtype: string containing the machine type
* mmodel: string containing the machine model
* mserial: string containing the machine serial number
* temp: array of ArcCheck measured temperatures, in C
* tilt: measured inclinometer tilt, in deg
* rot: measured inclinometer rotation, in deg
* dtype: string indicating diode type
* dosecal: absolute reference diode dose calibration, in Gy per count
* caltemp: temperature during calibration, in C
* tics: number of intervals/tics measured
* backgrounds: number of background measurements collected
* totaltime: total time between start/stop, in msec
* beamtime: total time where beam was measured, in msec
* z: array of detector spatial Z coordinates, in cm
* x: array of detector spatial X coordinates, in cm
* y: array of detector spatial Y coordinates, in cm
* theta: array of detector spatial cylindrical coordinate, in radians
* background: array of detector measured background, in counts/msec
* calibration: array of detector relative calibration
* data: 2D array of measured detector data

Below is an example of how this function is used:

```matlab
% Load SNC ArcCHECK data from two files
path = '/path/to/files/';
names = {'Head2_G90_to_G270_10deg.acm', 'Head3_G270_to_G90_10deg.acm'};
data = ParseSNCacm(path, names);

% Compute relative beam on fraction for measurement
fraction = data.beamtime / data.totaltime;
```

### ParseSNCprm

`ParseSNCprm()` extracts data from a SNC Profiler movie file (.prm) and returns the contents of the file as a MATLAB structure.  See below for a full list of the structure fields returned.  This function will display a progress bar while it loads (unless MATLAB was executed with the `-nodisplay`, `-nodesktop`, or `-noFigureWindows` flags).

The second input parameter can either contain one or multiple files.  If multiple files are selected, the returned data is a concatenation of all measured profiles.  However, only the header data (array calibration, background, ignore flags, etc) for the last file is returned.

If a given field is not found, this function will gracefully ingore it and the returned structure will not contain the field. If the field is found but no contents were specified, the returned field will be an empty cell array.

This function has been tested with SNC Profiler version 3.3.1 and .prm version 25. Note, there are additional fields in SNC .prm files that are not currently imported by this function. Additional fields can be added using the modular `search` variable, declared within this function. Refer to the documentation within the source code for more information.

The following variables are required for proper execution:

* path: string containing the path to the DICOM files
* names: string or cell array of strings containing the file(s) to be loaded

The following structure fields are returned upon successful completion:

* filerev: string containing the file revision letter
* filename: string containing the original filename
* timestamp: date and time .prm file was saved, as an integer
* description: string containing the SNC description
* institution: string containing the institution
* dcal: string containing the calibration file name
* version: string containing the SNC software version
* dorientation: string containing the detector orientation
* dssd: SSD, in cm
* dmodel: string containing the detector model
* dserial: string containing the detector serial number
* dfirmware: string containing the detector firmware
* dgain: gain
* dmode: string containing the measurement mode
* dinterval: collection interval, in ms
* mroom: string containing the room
* mtype: string containing the machine type
* mmodel: cstring containing the machine model
* mserial: string containing the machine S/N
* mbeamtype: string containing the beam type
* collimator: array of collimator values (left, right, top, bottom) in cm
* wangle: wedge angle
* mrate: array containing the dose rate, in MU/min
* mdose: array containing the dose delivered, in MU
* mangle: array containing the gantry angle
* cangle: array containing the collimator angle
* caltemp: calibration temperature, in C
* bpmpbt: array of Bp, Mp, and Bt calibration values
* dosecal: absolute reference chamber dose calibration, in Gy per count
* gain0: array of gain ratios for amp 0
* gain1: array of gain ratios for amp 1
* gain2: array of gain ratios for amp 2
* gain3: array of gain ratios for amp 3
* gain4: array of gain ratios for amp 4
* gain5: array of gain ratios for amp 5
* gain6: array of gain ratios for amp 6
* us1: array (dist, mv) of ultrasonic point 1 calibration
* us2: array (dist, mv) of ultrasonic point 2 calibration
* totaltime: total time between start/stop, in msec
* pulsestats: array of pulse statistics (idle, mmt, idle dur, mmt dur)
* num: array of number of detectors (X, Y, PD, ND, Ref, Z)
* detspacing: detector spacing, in cm
* background: array of detector measured background, in counts/timetic
* calibration: array of detector relative calibration
* ignore: array of detector ignore flags
* data: 2D array of measured detector data

Below is an example of how this function is used:

```matlab
% Load SNC Profiler data from one file
path = '/path/to/files/';
names = 'Head1_G90_27p3.prm';
data = ParseSNCprm(path, names);
```

### ParseSNCatxt

`ParseSNCtxt()` extracts data from a SNC Profiler ASCII File Export text file and returns the data returned as a MATLAB structure. See below for a full list of the structure fields returned.  This function will display a progress bar while it loads (unless MATLAB was executed with the `-nodisplay`, `-nodesktop`, or `-noFigureWindows` flags). Where indicated, returned arrays have a length equal to the number of profiles in the text file.  If a given field is not found, this function will gracefully ingore it and the returned structure will not contain the field. If the field is found but no contents were specified, the returned field will be an empty cell array.

This function has been tested with SNC Profiler version 3.3.1. Note, there are additional fields in SNC ASCII text files that are not currently imported by this function.  Additional fields can be added using the modular `search` variable, declared within this function. Refer to the documentation within the source code for more information.

The following variables are required for proper execution:

* path: string containing the path to the DICOM files
* name: string containing the file to be loaded

The following structure fields are returned upon successful completion:

* filenames: cell array of strings containing the filenames loaded
* timestamp: array of date and time that the file was saved, as integers
* description: cell array of strings containing the description
* institution: cell array of strings containing the institution
* version: cell array of strings containing the software version
* mroom: cell array of strings containing the room
* mtype: cell array of strings containing the machine type
* mmodel: cell array of strings containing the machine model
* mserial: cell array of strings containing the machine S/N
* mbeamtype: cell array of strings containing the beam type
* menergy: cell array of strings containing the beam energy
* wangle: cell array containing the wedge angle
* wtype: cell array of strings containing the wedge type
* mangle: array containing the gantry angle
* cangle: array containing the collimator angle
* cleft: array of the collimator left values, in cm
* cright: array of the collimator right values, in cm
* ctop: array of the collimator top values, in cm
* cbottom: array of the collimator bottom values, in cm
* mrate: array containing the dose rate, in MU/min
* mdose: array containing the dose delivered, in MU
* dorientation: cell array of strings containing the detector orientation
* dssd: array containing the SSD, in cm
* dcal: cell array of strings containing the calibration file
* dmodel: cell array of strings containing the detector model
* dserial: cell array of strings containing the detector serial number
* drev: cell array of strings containing the detector revision
* dfirmware: cell array of strings containing the detector firmware
* dmode: cell array of strings containing the measurement mode
* dgain: array containing the gain
* dinterval: array containing the collection interval, in ms
* cax: array of CAX/Normalized dose values, in cGy
* tics: array of timer tics, in microseconds
* xdata: 2D array of X axis data, where column one is the position (in cm), and columns 2:n+1 are the data for each measurement
* ydata: 2D array of Y axis data, where column one is the position (in cm), and columns 2:n+1 are the data for each measurement
* pdiag: 2D array of positive diagonal data, where column one is the position (in cm), and columns 2:n+1 are the data for each measurement
* ndiag: 2D array of negative diagonal data, where column one is the position (in cm), and columns 2:n+1 are the data for each measurement

Below is an example of how this function is used:

```matlab
% Load SNC ASCII data
path = '/path/to/files/';
name = 'Head1_G0.txt';
data = ParseSNCtxt(path, name);

% Plot Y axis profiles
figure;
hold on;
for i = 2:length(data.ydata)
plot(data.ydata{1}, data.ydata{i} * data.cax(i-1));
end
hold off;
```

## AnalyzeProfilerFields

`AnalyzeProfilerFields()` reads in a data structure (created via `ParseSNCprm()` or `ParseSNCtxt()`) and computes statistics for each field.  If multiple frames are detected when analyzing time-dependent ParseSNCprm data (indicated by gaps in time between measurements), each frame will be separated out and statistics will be computed independently.  The various statistics are returned as a MATLAB structure.

This function may be executed with a second input argument containing one or more reference profiles for each axis (xdata, ydata, pdiag, ndiag). If provided, additional statistics (such as Gamma comparisons) will be computed. If only one reference profile is provided, but multiple measured profiles are provided, the same reference profile will be applied to each measured profile. If multiple reference profiles are provided, the closest matching reference profile to the measured frame will be selected (using a correlation coefficient).

If reference data is provided, a second structure can also be returned by this function containing the same fields (and array sizes) as the first results array, but containing the statistics for the reference profile.

Finally, when correcting for ignored detectors or computing profile differences or gamma (which requires uniformly spaced data), interpolation may be required.  In these instances, spline-based interpolation is used.

This function has been tested with SNC Profiler version 3.3.1 and .prm version 25. 

The following variables are required for proper execution:

* varargin{1}: structure returned either by `ParseSNCprm()` or `ParseSNCtxt()` (see the documentation for these functions for more information on the fields contained)
* varargin{2} (optional): structure containing reference profile data. The following structure fields can be incuded: `xdata` (2 x n+1 array), `ydata` (2 x n+1 array), `pdiag` (2 x n+1 array), and `ndiag` (2 x n+1 array). To compute Gamma, must also include the fields `abs` (absolute gamma criterion), `dta` (distance to agreement, in cm), and `local` (flag indicating whether to perform local or global Gamma comparison).
* varargin{3} (optional): string indicating whether or not to normalize profiles when processing. Can be 'none', 'max', or 'center'.  If omitted, 'none' is assumed.  If reference data is not passed, this may also be passed as `varargin{2}` (see examples).

The following structure fields are returned for `varargout{1}` and `varargout{2}` upon successful completion:

* xdata: 2D array of background/ignore flag/calibration corrected X axis data, where column one is the position (in cm), and columns 2:n+1 are the dose for each frame/profile (in Gy)
* ydata: 2D array of background/ignore flag/calibration corrected Y axis data, where column one is the position (in cm), and columns 2:n+1 are the dose for each frame/profile (in Gy)
* pdiag: 2D array of background/ignore flag/calibration corrected positive diagonal axis data, where column one is the position (in cm), and columns 2:n+1 are the dose for each frame/profile (in Gy)
* ndiag: 2D array of background/ignore flag/calibration corrected negative diagonal axis data, where column one is the position (in cm), and columns 2:n+1 are the dose for each frame/profile (in Gy)
* tdata (optional): 2D array of time-dependent detector channel data, where column one is the absolute time and column two is the differential central channel response. If time dependent data is not available (i.e. ParseSNCtxt data), this field is not returned.
* xfwhm: vector of X axis Full Width at Half Maximum(s) for each field. Note, if profile edges cannot be found for a given profile, the FWHM and edges values will be zero.
* xedges: n x 2 array of left and right FWHM-defined X axis field edges
* yfwhm: vector of Y axis Full Width at Half Maximum(s) for each field
* yedges: n x 2 array of left and right FWHM-defined Y axis field edges
* pfwhm: vector of positive diagonal axis Full Width at Half Maximum(s) for each field
* pedges: n x 2 array of left and right FWHM-defined positive diagonal axis field edges
* nfwhm: vector of negative diagonal axis Full Width at Half Maximum(s) for each field
* nedges: n x 2 array of left and right FWHM-defined negative diagonal axis field edges
* xflat: vector of central 80% flatness for X axis profiles
* xsym: vector of areal symmetry for X axis profiles
* yflat: vector of central 80% flatness for Y axis profiles
* ysym: vector of areal symmetry for Y axis profiles
* pflat: vector of central 80% flatness for positive diagonal profiles
* psym: vector of areal symmetry for positive diagonal profiles
* nflat: vector of central 80% flatness for negative diagonal profiles
* nsym: vector of areal symmetry for negative diagonal profiles

The following additional structure fields are returned in `varargout{1}` if reference data was provided and `nargout == 2`:

* ref (optional): vector indicating which reference profile was selected
* corr (optional): 4 x n x m 3D array containing the correlation coefficient between each of n measured profiles, m reference profiles, and axis (x, y, p, n)
* xdiff (optional): 2D array of X axis differences, where column one is the position (in cm), and columns 2:n+1 are the abs differences
* ydiff (optional): 2D array of Y axis differences, where column one is the position (in cm), and columns 2:n+1 are the abs differences
* pdiff (optional): 2D array of positive diagonal differences, where column one is the position (in cm), and columns 2:n+1 are the abs differences
* ndiff (optional): 2D array of negative diagonal differences, where column one is the position (in cm), and columns 2:n+1 are the abs differences
* xgamma (optional): 2D array of X axis Gamma values, where column one is the position (in cm), and columns 2:n+1 are the Gamma indices
* ygamma (optional): 2D array of Y axis Gamma values, where column one is the position (in cm), and columns 2:n+1 are the Gamma indices
* pgamma (optional): 2D array of positive diagonal Gamma values, where column one is the position (in cm), and columns 2:n+1 are the Gamma indices
* ngamma (optional): 2D array of negative diagonal Gamma values, where column one is the position (in cm), and columns 2:n+1 are the Gamma indices

Below is an example of how this function is used:

```matlab
% Load SNC Profiler data from one file
path = '/path/to/files/';
name = 'Head1_27p3.prm';
data = ParseSNCprm(path, name);

% Compute statistics on Profiler data (without normalization)
results = AnalyzeProfilerFields(data, 'none');

% Load reference profiles from a DICOM dose file
refdata = LoadProfilerDICOMReference(...
'AP_27P3X27P3_PlaneDose_Vertical_Isocenter.dcm');

% Set Gamma criteria
refdata.abs = 2;    % percent
refdata.dta = 0.1;  % cm
refdata.local = 0;  % perform global analysis

% Compute statistics again, this time comparing to reference criteria
% and normalizing profiles to Profiler center
[results, refresults] = AnalyzeProfilerFields(data, refdata, 'center');

% Plot X axis Gamma index
figure;
subplot(2,2,1);
hold on;
plot(results.xdata(1,:), results.xdata(2,:));
plot(refresults.xdata(1,:), refresults.xdata(2,:));
plot(results.xgamma(1,:), results.xgamma(2,:));
hold off;
title(sprintf('xgamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
xlabel('x axis (cm)');

% Plot Y axis Gamma index
subplot(2,2,2);
hold on;
plot(results.ydata(1,:), results.ydata(2,:));
plot(refresults.ydata(1,:), refresults.ydata(2,:));
plot(results.ygamma(1,:), results.ygamma(2,:));
hold off;
title(sprintf('ygamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
xlabel('y axis (cm)');

% Plot positive diagonal axis Gamma index
subplot(2,2,3);
hold on;
plot(results.pdiag(1,:), results.pdiag(2,:));
plot(refresults.pdiag(1,:), refresults.pdiag(2,:));
plot(results.pgamma(1,:), results.pgamma(2,:));
hold off;
title(sprintf('pgamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
xlabel('positive diagonal (cm)');

% Plot negative diagonal axis Gamma index
subplot(2,2,4);
hold on;
plot(results.ndiag(1,:), results.ndiag(2,:));
plot(refresults.ndiag(1,:), refresults.ndiag(2,:));
plot(results.ngamma(1,:), results.ngamma(2,:));
hold off;
title(sprintf('ngamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
xlabel('negative diagonal (cm)');
```

## AnalyzeACFields

`AnalyzeACFields()` reads in an ArcCHECK data structure (created via `ParseSNCacm()`), identifies individual fields, and computes the FWHM-defined center of each field (in cylindrical coordinates).   Each field detector data is interpolated to a 2D cylindrical frame. The results are returned as a MATLAB structure.

This function will display a progress bar while it loads (unless MATLAB was executed with the `-nodisplay`, `-nodesktop`, or `-noFigureWindows` flags).

This structure uses a cylindrical coordinate system to define ArcCHECK data, where Y is positioned along the central long axis of the ArcCHECK and theta is the angle (in degrees) from 0, where 0 is the ArcCHECK zero position (typically positioned in the positive IEC-Z direction)

This function has been tested with SNC Patient version 6.2.3 and .acm file rev C. 

When computing the 2D cylindrical frame, the MATLAB `scatteredInterpolant` class is used with a bilinear method.  When computing the FWHM-defined edges, splines-based interpolation is applied.

The following variables are required for proper execution:

* varargin{1}: structure returned either by `ParseSNCacm()` (see the documentation for this function for more information on the fields 
contained)
* varargin{2} (optional): number indicating the angle (in degrees) to offset the measured angles, to account for ArcCheck roll

The following structure fields are returned upon successful completion:

* detectors: 1386 x n array of detector data, corrected for background, relative calibration, and absolute calibration (in Gy)
* frames: 360 x 201 x n array of n fields, where each 2D array is an
* interpolated cylindrical map (theta x Y) of measured dose (in Gy)
* alpha: 2 x n array of entrance and exit FWHM-defined profile center theta values (in degrees). Note, if profile edges cannot be found for a given profile, the alpha and beta values will be zero.
* beta: 2 x n array of entrance and exit FWHM-defined profile center Y axis values (in cm)

Below is an example of how this function is used:

```matlab
% Load SNC ArcCHECK data from one file
path = '/path/to/files/';
names = 'Head3_G270_to_G90_10deg.acm';
data = ParseSNCacm(path, names);

% Analyze ArcCHECK data
results = AnalyzeACFields(data);

% Loop through frames, plotting result
figure;
for i = 1:size(results.frames, 3)

    % Plot frame
    imagesc(circshift(results.frames(:,:,i), -180, 2));
    set(gca,'XTick', 1:30:361);
    set(gca,'XTickLabel', -180:30:180);
    xlabel('ArcCHECK Angle (deg)');
    set(gca,'YTick', 1:20:201);
    set(gca,'YTickLabel', 10:-2:-10);
    ylabel('ArcCHECK Y (cm)');
    title(sprintf('Frame %i', i));

    % Update plot and pause temporarily
    drawnow;
    pause(0.1);
end
```

## Event Calling

These functions optionally return execution status and error information to an `Event()` function. If available in the MATLAB path, `Event()` will be called with one or two variables: the first variable is a string containing the status information, while the second is the status classification (WARN or ERROR). If the status information is only informative, the second argument is not included.  Finally, if no `Event()` function is available errors will still be thrown via the standard `error()` MATLAB function.

## Gamma Computation Methods

When comparing measured to reference profiles, a Gamma analysis is performed based on the formalism presented by D. A. Low et. al., [A technique for the quantitative evaluation of dose distributions.](http://www.ncbi.nlm.nih.gov/pubmed/9608475), Med Phys. 1998 May; 25(5): 656-61.  In this formalism, the Gamma quality index *&gamma;* is defined as follows for each point in measured dose/response volume *Rm* given the reference dose/response volume *Rc*:

*&gamma; = min{&Gamma;(Rm,Rc}&forall;{Rc}*

where:

*&Gamma; = &radic; (r^2(Rm,Rc)/&Delta;dM^2 + &delta;^2(Rm,Rc)/&Delta;DM^2)*,

*r(Rm,Rc) = | Rc - Rm |*,

*&delta;(Rm,Rc) = Dc(Rc) - Dm(Rm)*,

*Dc(Rc)* and *Dm(Rm)* represent the reference and measured doses at each *Rc* and *Rm*, respectively, and

*&Delta;dM* and *&Delta;DM* represent the absolute and Distance To Agreement Gamma criterion (by default 3%/3mm), respectively.  

The absolute criterion is typically given in percent and can refer to a percent of the maximum dose (commonly called the global method) or a percentage of the voxel *Rm* being evaluated (commonly called the local method).  The application is capable of computing gamma using either approach, and can be set in `AnalyzeProfilerFields()` by passing the input argument `varargin{2}.local` as 0 or 1.  If not provided, the global method (0) is applied.

To improve computation efficiency, the computation space *&forall;{Rc}* is limited to twice the distance to agreement parameter.  Thus, the maximum "real" Gamma index returned by the application is 2.
