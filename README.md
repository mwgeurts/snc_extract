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
* [Event Calling](README.md#event-calling)
 
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

This function has been tested with SNC Profiler version 3.3.1 and .prm version 25. Note, there are additional fields in SNC .prm files that are not currently imported by this function. Additional fields can be added using the modular search variable, declared within this function. Refer to the documentation within the source code for more information.

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
* numdetectors: array of number of detectors (X, Y, PD, ND, Ref, Z)
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

This function has been tested with SNC Profiler version 3.3.1. Note, there are additional fields in SNC ASCII text files that are not currently imported by this function.  Additional fields can be added using the modular search variable, declared within this function. Refer to the documentation within the source code for more information.

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

## Event Calling

These functions optionally return execution status and error information to an `Event()` function. If available in the MATLAB path, `Event()` will be called with one or two variables: the first variable is a string containing the status information, while the second is the status classification (WARN or ERROR). If the status information is only informative, the second argument is not included.  Finally, if no `Event()` function is available errors will still be thrown via the standard `error()` MATLAB function.
