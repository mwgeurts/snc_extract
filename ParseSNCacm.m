function data = ParseSNCacm(path, names)
% ParseSNCacm extracts data from a SNC ArcCHECK movie file (.acm) and
% returns the contents of the file as a MATLAB structure.  See below for a 
% full list of the structure fields returned.  This function will display a 
% progress bar while it loads (unless MATLAB was executed with the 
% -nodisplay, -nodesktop, or -noFigureWindows flags).
%
% The second input parameter can either contain one or multiple files.  If
% multiple files are selected, the returned data is a concatenation of all
% measured profiles.  However, only the header data (array calibration, 
% background, ignore flags, etc) for the last file is returned.
%
% If a given field is not found, this function will gracefully ingore it 
% and the returned structure will not contain the field. If the
% field is found but no contents were specified, the returned field will be
% an empty cell array.
%
% This function has been tested with SNC Patient version 6.2.3 and .acm 
% file rev C. Note, there are additional fields in SNC .acm files that are 
% not currently imported by this function. Additional fields can be added 
% using the modular search variable, declared within this function. Refer 
% to the documentation within the source code for more information.
%
% The following variables are required for proper execution:
%   path: string containing the path to the DICOM files
%   names: string or cell array of strings containing the file(s) to be 
%       loaded
%
% The following structure fields are returned upon successful completion:
%   filerev: string containing the file revision letter
%   filename: string containing the original filename
% 	timestamp: date and time .acm file was saved, as an integer
%   institution: string containing the institution
%   dcal: string containing the calibration file name
%   version: string containing the SNC software version
%   dorientation: string containing the ArcCHECK orientation
%   ssd: SSD, in cm
%   shiftx: phantom X shift, in mm
%   shifty: phantom Y shift, in mm
%   shiftz: phantom Z shift, in mm
%   rotx: phantom X rotation, in deg
%   roty: phantom Y rotation, in deg
%   rotz: phantom Z rotation, in deg
%   dmodel: string containing the detector model
%   dserial: string containing the detector serial number
%   dfirmware: string containing the detector firmware revision number
%   dinterval: collection interval, in ms
%   dthreshold: background threshold, in percent
%   dplug: string indicating whether plug was present for acquisition
%   dplugdose: measured dose (in Gy) to central cavity, if provided
%   mroom: string containing the treatment room
%   mtype: string containing the machine type
%   mmodel: string containing the machine model
%   mserial: string containing the machine serial number
%   temp: array of ArcCheck measured temperatures, in C
%   tilt: measured inclinometer tilt, in deg
%   rot: measured inclinometer rotation, in deg
%   dtype: string indicating diode type
%   dosecal: absolute reference diode dose calibration, in Gy per count
%   caltemp: temperature during calibration, in C
%   tics: number of intervals/tics measured
%   backgrounds: number of background measurements collected
%   totaltime: total time between start/stop, in msec
%   beamtime: total time where beam was measured, in msec
%   z: array of detector spatial Z coordinates, in cm
%   x: array of detector spatial X coordinates, in cm
%   y: array of detector spatial Y coordinates, in cm
%   theta: array of detector spatial cylindrical coordinate, in radians
%   background: array of detector measured background, in counts/msec
%   calibration: array of detector relative calibration
%   data: 2D array of measured detector data
%
% Below is an example of how this function is used:
%
%   % Load SNC ArcCHECK data from two files
%   path = '/path/to/files/';
%   names = {'Head2_G90_to_G270_10deg', 'Head3_G270_to_G90_10deg'};
%   data = ParseSNCacm(path, names);
%
%   % Compute relative beam on fraction for measurement
%   fraction = data.beamtime / data.totaltime;
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

% If not cell array, cast as one
if ~iscell(names); names = cell({names}); end

% Log start of file load and start timer
if exist('Event', 'file') == 2
    Event(['Loading ACM file ', strjoin(names, '\nLoading ACM file ')]);
    tic;
end

% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    
    % Start waitbar
    progress = waitbar(0, 'Parsing SNC ArcCHECK file(s)');
end

% Declare search variables. This array specifies what lines are extracted
% from the file, and into what format. The first column is the stored 
% variable name, the second is the search string, and the third is the data 
% type. Currently supported type values are string, float, datenum, csv, 
% position, vector, and data; see while loop for specifics.
search = {
    'filerev'   'File Revision:'    'string'
    'filename'  'Filename:' 'string'
    'timestamp' 'Date:' 'datenum'
    'institution'   'Institution:'  'string'
    'dcal'   'Calibration File:' 'string'
    'version'   'Software Version:' 'string'
    'dorientation'   'Orientation:'  'string'
    'ssd'   'SSD:'  'float'
    'shiftx'    'Shift X:'  'float'
    'shifty'    'Shift Y:'  'float'
    'shiftz'    'Shift Z:'  'float'
    'rotx'  'Rotation X:' 'float'
    'roty'  'Rotation Y:' 'float'
    'rotz'  'Rotation Z:' 'float'
    'dmodel' 'Collector Model:' 'string'
    'dserial'   'Collector Serial:' 'string'
    'dfirmware' 'Firmware Version:' 'string'
    'dinterval' 'Collection Interval:'  'float'
    'dthreshold'    'Background Threshold:' 'float'
    'dplug' 'Plug Present:' 'string'
    'dplugdose' 'Measured Cavity Dose:' 'float'
    'mroom' 'Room:' 'string'
    'mtype' 'Machine Type:' 'string'
    'mmodel'    'Machine Model:'    'string'
    'mserial'   'Machine Serial Number:'    'string'
    'temp'  'Temperature:'  'csv'
    'tilt'  'Inclinometer Tilt:' 'float'
    'rot'   'Inclinometer Rotation:'    'float'
    'dtype' 'Diode Type:'   'string'
    'dosecal'   'Dose Per Count:'   'float'
    'caltemp'   'Temperature During Calibration:'   'float'
    'tics'  'Number of Updates:'    'float'
    'backgrounds'   'Number of Backgrounds:'    'float'
    'totaltime'  'Total Time:'   'float'
    'beamtime'  'Beam ON Time:' 'float'
    'z' 'Detector Spacing:' 'vector'
    'x' 'Concatenation:'    'vector'
    'y' 'Imported Data:'    'vector'
    'background'    'Background'    'vector'
    'calibration'   'Calibration'   'vector'  
    'data'  'Data:' 'data'
};

% Initialize return variable
data = struct();

% Execute in try/catch statement
try

% Loop through each file, concatenating data arrays
for i = 1:length(names)
    
    % Attempt to open file handle to data
    fid = fopen(fullfile(path, names{i}), 'r');

    % Verify file handle is valid
    if fid >= 3
        if exist('Event', 'file') == 2
            Event('Read handle successfully established');
        end
    else
        if exist('Event', 'file') == 2
            Event(['Read handle not successful for ', names{i}], 'ERROR');
        else
            error(['Read handle not successful for ', names{i}]);
        end
    end
    
    % Retrieve the first line in the file
    tline = fgetl(fid);

    % If file does not start with 'ArcCHECK Measurement file'
    if ~strcmp(tline(1:25), 'ArcCHECK Measurement file')

        % The file may not be in correct format, so throw an error
        if exist('Event', 'file') == 2
            Event('File is not in expected format', 'ERROR');
        else
            error('File is not in expected format');
        end
    end

    % While the end-of-file has not been reached
    while ~feof(fid)
    
        % Retrieve the next line in the file
        tline = fgetl(fid);

        % Loop through each search variable
        for j = 1:size(search, 1)
            
            % If search variable is found
            if length(tline) >= length(char(search(j,2)))+1 && ...
                    strcmp(sprintf('%s\t', char(search(j,2))), ...
                    tline(1:length(char(search(j,2)))+1))
                
                % Update waitbar
                if exist('progress', 'var') && ishandle(progress)
                    waitbar(((i-1) * size(search, 1) + j - 1) / ...
                        (size(search, 1) * length(names)), progress);
                end
                
                % If returning a string
                if strcmp(search(j,3), 'string')
                    
                    % Store results as a string
                    data.(char(search(j,1))) = char(regexp(tline(length(...
                        char(search(j,2)))+2:end), '^([^\t]+)', 'match'));
                
                % Otherwise, if returning a float
                elseif strcmp(search(j,3), 'float')
                
                    % Store results as a double
                    data.(char(search(j,1))) = str2double(regexp(tline(...
                        length(char(search(j,2)))+2:end), '^([^\t]+)', ...
                        'match'));
                
                % Otherwise, if returning a datenum
                elseif strcmp(search(j,3), 'datenum')
                
                    % Temporarily store cell array
                    C = regexp(tline(length(char(search(j,2)))+2:end), ...
                        '^([^\t]+)\tTime:\t([^\t]+)', 'tokens');
                
                    % Compute and store datenum
                    data.(char(search(j,1))) = datenum(strjoin(C{1}, ' '));
                    
                    % Clear temporary variable
                    clear C;
                 
                % Otherwise, if returning a csv array
                elseif strcmp(search(j,3), 'csv')
                    
                    % Store results as a double
                    data.(char(search(j,1))) = str2double(strsplit(...
                        char(regexp(tline(length(char(search(j,2)))+2:end), ...
                        '^([^\t]+)', 'match')), ','));
                    
                % Otherwise, if returning a vector
                elseif strcmp(search(j,3), 'vector')
                    
                    % Temporarily store cell array
                    C = str2double(strsplit(tline, '\t'));
                    
                    % Store values, removing NaNs
                    data.(char(search(j,1))) = single(C(~isnan(C)));
                    
                    % Clear temporary variable
                    clear C;
                    
                % Otherwise, if returning a data array
                elseif strcmp(search(j,3), 'data')
                    
                    % If return structure field does not already exist,
                    % initialize it (for concatentation)
                    if ~isfield(data, char(search(j,1)))
                        data.(char(search(j,1))) = [];
                    end
                    
                    % Determine number of data elements
                    n = length(strsplit(tline, '\t')) - 1;
                    
                    % Move file pointer back to beginning of line
                    if fseek(fid, -length(tline), 0) == 0
                    
                        % Textscan remaining lines in file, storing results
                        % to temporary array (given number of elements
                        % determined above)
                        C = textscan(fid, ['%s', repmat(' %f', 1, n)]);
                        
                        % Remove string data from first column (to allow
                        % cell2mat)
                        C{1,1} = zeros(size(C{1,2},1),1);
                        
                        % Concatenate new array onto any existing data
                        data.(char(search(j,1))) = vertcat(...
                            data.(char(search(j,1))), single(cell2mat(C)));

                    % Otherwise throw an error
                    else
                        if exist('Event', 'file') == 2
                            Event('Error moving file pointer', 'ERROR');
                        else
                            error('Error moving file pointer');
                        end
                    end
                    
                    % Clear temporary variables
                    clear C n;
                    
                % Otherwise, return an error
                else
                    if exist('Event', 'file') == 2
                        Event('Search variable type is not supported', ...
                            'ERROR');
                    else
                        error('Search variable type is not supported');
                    end
                end
            end
        end
    end
    
    % Close file handle
    fclose(fid);
    
    % Clear temporary files
    clear fid tline j;
end

% Log SNC file version, if available
if exist('Event', 'file') == 2 && isfield(data, 'filerev') && ...
        ~isempty(data.filerev)
    Event(sprintf('SNC ACM File Revision %s', data.filerev));
end

% Log SNC application version, if available
if exist('Event', 'file') == 2 && isfield(data, 'version') && ...
        ~isempty(data.version)
    Event(sprintf('SNC Patient Version %s', data.version));
end

% If detector X and Z data was found
if isfield(data, 'z') && isfield(data, 'x')
    
    % Log action
    if exist('Event', 'file') == 2
        Event('Computing detector cylindrical coordinates');
    end
    
    % Compute cylindrical coordinates
    data.theta = atan2d(data.z, data.x) + 90;
end

% Log completion of function
if exist('Event', 'file') == 2
    Event(sprintf('Successfully parsed file(s) in %0.3f seconds', toc));
end

% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

% Clear temporary files
clear progress search i;

% Catch errors, log, and rethrow
catch err
    
    % Delete progress handle if it exists
    if exist('progress', 'var') && ishandle(progress), delete(progress); end
    
    % Log error
    if exist('Event', 'file') == 2
        Event(getReport(err, 'extended', 'hyperlinks', 'off'), 'ERROR');
    else
        rethrow(err);
    end
end