function data = LoadProfilerDICOMReference(varargin)
% LoadProfilerDICOMReference reads in a list of 2D DICOM RT Dose files 
% (file names with path) and extracts profiles along the locations of the
% SNC Profiler axes (X, Y, and positive/negative diagonals).  The profiles
% are returned as a structure that is compatible with AnalyzeProfilerFields
% (refer to the documentation in the README for more information on how
% reference data can be used by this function).
%
% This function uses MLATAB's dicominfo and dicomread functions to extract
% the file contents of input file.  This function assumes that the Profiler
% is aligned with the Y axis along the DICOM Y axis and the X axis along
% with DICOM X axis.
% 
% If the voxel sizes or start coordinates differ when loading multiple
% files, this function will resample all datasets to the resolution and 
% frame of reference of the first set provided. 
%
% The following variables are required for proper execution:
%   varargin: cell array of strings containing the full path and file name 
%       of each DICOM RT Dose file to be loaded.
%
% The following structure fields are returned upon successful completion:
%   xdata: 2D array of X axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the data for each input
%   ydata: 2D array of Y axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the data for each input
%   pdiag: 2D array of positive diagonal data, where column one is the
%       position (in cm), and columns 2:n+1 are the data for each input
%   ndiag: 2D array of negative diagonal data, where column one is the
%       position (in cm), and columns 2:n+1 are the data for each input
%
% Below is an example of how this function is used:
%
%   % Load two DICOM reference profiles
%   refdata = LoadProfilerDICOMReference(...
%       'AP_27P3X27P3_PlaneDose_Vertical_Isocenter.dcm', ...
%       'AP_10P5X10P5_PlaneDose_Vertical_Isocenter.dcm');
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

% Execute in try/catch statement
try
    
% Log start
if exist('Event', 'file') == 2
    Event('Loading SNC Profiler reference datasets');
    tic;
end

% Initialize return variable
data = struct;
    
% Loop through input arguments
for i = 1:nargin
    
    % Log file being loaded
    if exist('Event', 'file') == 2
        Event(sprintf('Loading %s', varargin{i}));
    end
    
    % Load DICOM data
    info = dicominfo(varargin{i});
    width = info.PixelSpacing/10;
    start = [info.ImagePositionPatient(3); info.ImagePositionPatient(1)]/10;
    image = single(dicomread(varargin{i})) * info.DoseGridScaling;
    
    % Generate mesh
    [meshX, meshY]  = meshgrid(start(2):width(2):start(2) + width(2) * ...
        (size(image, 2) - 1), start(1) + width(1) * ...
        (size(image, 1) - 1): -width(1):start(1));
    
    % If a frame of reference does not yet exist
    if ~isfield(data, 'xdata')
       
        % Log file being loaded
        if exist('Event', 'file') == 2
            Event('Generating return data frame of reference');
        end
        
        % Store sorted profile positions from resolution of first file
        data.xdata(1,:) = sort(meshX(1,:));
        data.ydata(1,:) = sort(meshY(:,1));
        
        % Store diagonal positions as the lerger of the x/y dimensions
        if size(data.xdata, 2) < size(data.ydata, 2)
            data.pdiag(1,:) = data.ydata(1,:) * sqrt(2);
            data.ndiag(1,:) = data.ydata(1,:) * sqrt(2);
        else
            data.pdiag(1,:) = data.xdata(1,:) * sqrt(2);
            data.ndiag(1,:) = data.xdata(1,:) * sqrt(2);
        end
        
    end
    
    % Log file being loaded
    if exist('Event', 'file') == 2
        Event('Interpolating SNC Profiler axes');
    end
    
    % Interpolate each axis and store reference profile
    data.xdata(i+1,:) = interp2(meshX, meshY, image, data.xdata(1,:), ...
        zeros(1, size(data.xdata, 2)), '*linear', 0);
    data.ydata(i+1,:) = interp2(meshX, meshY, image, zeros(1, ...
        size(data.ydata, 2)), data.ydata(1,:), '*linear', 0);
    data.pdiag(i+1,:) = interp2(meshX, meshY, image, data.pdiag(1,:) ...
        * sqrt(2)/2, data.pdiag(1,:) * sqrt(2)/2, '*linear', 0);
    data.ndiag(i+1,:) = interp2(meshX, meshY, image, data.ndiag(1,:) ...
        * sqrt(2)/2, -data.ndiag(1,:) * sqrt(2)/2, '*linear', 0);
    
    % Clear temporary variables
    clear info width start image meshX meshY;
end

% Clear temporary variables
clear i;

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Reference datasets successfully loaded in ', ...
        '%0.3f seconds'], toc));
end

% Catch errors, log, and rethrow
catch err
    if exist('Event', 'file') == 2
        Event(getReport(err, 'extended', 'hyperlinks', 'off'), 'ERROR');
    else
        rethrow(err);
    end
end


