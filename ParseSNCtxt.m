function data = ParseSNCtxt(path, name)
% ParseSNCtxt extracts data from a SNC Profiler ASCII File Export text file 
% and returns the data returned as a MATLAB structure. See below for a full
% list of the structure fields returned.  This function will display a 
% progress bar while it loads (unless MATLAB was executed with the 
% -nodisplay, -nodesktop, or -noFigureWindows flags). Where indicated, 
% returned arrays have a length equal to the number of profiles in the text 
% file.  If a given field is not found, this function will gracefully
% ingore it and the returned structure will not contain the field.
%
% This function has been tested with SNC Profiler version 3.3.1.
%
% The following variables are required for proper execution:
%   path: string containing the path to the DICOM files
%   name: string containing the file to be loaded
%
% The following structure fields are returned upon successful completion:
%   timestamp: array of date and time that the file was saved, as integers
%   description: array of strings containing the description
%   institution: array of strings containing the institution
%   version: array of strings containing the software version
%   mtype: array of strings containing the machine type
%   mmodel: array of strings containing the machine model
%   mserial: array of strings containing the machine serial number
%   beamtype: array of strings containing the beam type
%   energy: array of strings containing the beam energy
%   wangle: array containing the wedge angle
%   wtype: array of strings containing the wedge type
%   cangle: array of strings containing the collimator angle
%   cleft: array containing the collimator left values, in cm
%   cright: array containing the collimator right values, in cm
%   ctop: array containing the collimator top values, in cm
%   cbottom: array containing the collimator bottom values, in cm
%   rate: array containing the dose rate, in MU/min
%   dose: array containing the dose delivered, in MU
%   orientation: array of strings containing the detector orientation
%   ssd: array containing the SSD, in cm
%   cal: array of strings containing the calibration file used
%   dmodel: array of strings containing the detector model
%   dserial: array of strings containing the detector serial number
%   drev: array of strings containing the detector revision
%   dfirm: array of strings containing the detector firmware
%   mode: array of strings containing the measurement mode
%   gain: array containing the gain
%   interval: array containing the collection interval, in ms
%   tics: array of timer tics
%   xdata: 2-D array of X axis positions, where column one is the position
%       (in cm), and columns 2:n+1 are the data for each measurement
%   ydata: 2-D array of Y axis positions, where column one is the position
%       (in cm), and columns 2:n+1 are the data for each measurement
%   pdiag: 2-D array of positive diagonal positions, where column one is 
%       the position (in cm), and columns 2:n+1 are the data for each 
%       measurement
%   ndiag: 2-D array of negative diagonal positions, where column one is 
%       the position (in cm), and columns 2:n+1 are the data for each 
%       measurement
%
% Below is an example of how this function is used:
%
%   % Load SNC ASCII data
%   path = '/path/to/files/';
%   name = 'Head1_G0.txt';
%   data = ParseSNCtxt(path, name);
%
%   % Plot X axis profiles
%   figure;
%   hold on;
%   for i = 2:length(data.xdata)
%       plot(data.xdata(1,:), data.xdata(1,:));
%   end
%   hold off;
%
% Copyright (C) 2015 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.
