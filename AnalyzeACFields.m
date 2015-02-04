function results = AnalyzeACFields(varargin)
% AnalyzeACFields reads in an ArcCHECK data structure (created via 
% ParseSNCacm), identifies individual fields, and computes the FWHM-defined
% center of each field (in cylindrical coordinates).   Each field detector 
% data is interpolated to a 2D cylindrical frame. The results are returned 
% as a MATLAB structure.
%
% This structure uses a cylindrical coordinate system to define ArcCHECK 
% data, where Y is positioned along the central long axis of the ArcCHECK
% and theta is the angle (in degrees) from 0, where 0 is the ArcCHECK zero
% position (typically positioned in the positive IEC-Z direction)
%
% This function has been tested with SNC Patient version 6.2.3 and .acm 
% file rev C. 
%
% The following variables are required for proper execution:
%   varargin{1}: structure returned either by ParseSNCacm (see the 
%       documentation for this function for more information on the fields 
%       contained)
%   varargin{2} (optional): number indicating the angle (in degrees) to 
%       offset the measured angles, to account for ArcCheck roll
%
% The following structure fields are returned upon successful completion:
%   frames: 360 x 201 x n array of n fields, where each 2D array is an
%   	interpolated cylindrical map (theta x Y) of measured dose
%   alpha: 2 x n array of entrance and exit FWHM-defined profile center
%       theta values (in degrees)
%   beta: 2 x n array of entrance and exit FWHM-defined profile center Y 
%       axis values (in cm)
%
% Below is an example of how this function is used:
%
%
%
%
%
%
%
%
%
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

%% Validate inputs
% Check the number of inputs
if nargin == 0 || nargin > 2
    if exist('Event', 'file') == 2
        Event('Incorrect number of arguments passed to function', 'ERROR');
    else
        error('Incorrect number of arguments passed to function');
    end
end

% If first argument contains the correct fields
if isstruct(varargin{1}) && isfield(varargin{1}, 'data') && ...
        isfield(varargin{1}, 'y') && ...
        isfield(varargin{1}, 'theta') && ...
        isfield(varargin{1}, 'background') && ...
        isfield(varargin{1}, 'calibration')
    
% Otherwise, type is invalid so throw an error
else
    if exist('Event', 'file') == 2
        Event('Invalid first input argument', 'ERROR');
    else
        error('Invalid first input argument');
    end
end

% Check for rotation
if nargin == 2 && isnumeric(varargin{2})
    
    % Set rotation variable
    rot = varargin{2};
    
    % Log event
    if exist('Event', 'file') == 2
        Event(sprintf('ArcCHECK angles adjusted by %g deg', rot));
    end
    
% Otherwise do not rotate profiles
else
    rot = 0;
end

% Execute in try/catch statement
try

% Log start of analysis and start timer
if exist('Event', 'file') == 2
    Event('Analyzing SNC ArcCHECK data');
    timer = tic;
end 

% Calculate interpolation meshgrids
[itheta, iY] = meshgrid(0:359, -10:0.1:10);

% Log interpolation grid settings
if exist('Event', 'file') == 2
    Event('Interpolation grid set to [0:359 deg, -10:0.1:10 cm]');
end



















% Log completion
if exist('Event', 'file') == 2
    
    % Log event 
    Event(sprintf(['SNC ArcCHECK analysis successfully completed in ', ...
      '%0.3f seconds'], toc(timer)));
  
    % Clear temporary variable
    clear timer;
end

% Catch errors, log, and rethrow
catch err
    if exist('Event', 'file') == 2
        Event(getReport(err, 'extended', 'hyperlinks', 'off'), 'ERROR');
    else
        rethrow(err);
    end
end
