function data = ParseSNCtxt(path, name, varargin)
% ParseSNCtxt extracts data from a SNC Profiler ASCII File Export text file 
% and returns the data returned as a MATLAB structure. See below for a full
% list of the structure fields returned.  This function will display a 
% progress bar while it loads (unless MATLAB was executed with the 
% -nodisplay, -nodesktop, or -noFigureWindows flags). Where indicated, 
% returned arrays have a length equal to the number of profiles in the text 
% file.  If a given field is not found, this function will gracefully
% ingore it and the returned structure will not contain the field. If the
% field is found but no contents were specified, the returned field will be
% an empty cell array.
%
% This function has been tested with SNC Profiler version 3.3.1. Note,
% there are additional fields in SNC ASCII text files that are not
% currently imported by this function.  Additional fields can be added
% using the modular search variable, declared within this function. Refer
% to the documentation within the source code for more information.
%
% The following variables are required for proper execution:
%   path: string containing the path to the TXT files
%   name: string containing the file to be loaded
%   varargin: optional name/value pairs 'Type' and 'Progress'. If Type is 
%       not specifiedthe function will attempt to recognize it). Type can
%       be 'PROFILER' or 'ARCCHECK'. Progress is a boolean, and will
%       determine whether or not to display the progress bar.
%
% The following structure fields are returned upon successful completion:
%   filenames: cell array of strings containing the filenames loaded
%   timestamp: array of date and time that the file was saved, as integers
%   description: cell array of strings containing the description
%   institution: cell array of strings containing the institution
%   version: cell array of strings containing the software version
%   mroom: cell array of strings containing the room
%   mtype: cell array of strings containing the machine type
%   mmodel: cell array of strings containing the machine model
%   mserial: cell array of strings containing the machine S/N
%   mbeamtype: cell array of strings containing the beam type
%   menergy: cell array of strings containing the beam energy
%   wangle: cell array containing the wedge angle
%   wtype: cell array of strings containing the wedge type
%   mangle: array containing the gantry angle
%   cangle: array containing the collimator angle
%   cleft: array of the collimator left values, in cm
%   cright: array of the collimator right values, in cm
%   ctop: array of the collimator top values, in cm
%   cbottom: array of the collimator bottom values, in cm
%   mrate: array containing the dose rate, in MU/min
%   mdose: array containing the dose delivered, in MU
%   dorientation: cell array of strings containing the detector 
%       orientation
%   dssd: array containing the SSD, in cm
%   dcal: cell array of strings containing the calibration file
%   dmodel: cell array of strings containing the detector model
%   dserial: cell array of strings containing the detector serial 
%       number
%   drev: cell array of strings containing the detector revision
%   dfirmware: cell array of strings containing the detector 
%       firmware
%   dmode: cell array of strings containing the measurement mode
%   dgain: array containing the gain
%   dinterval: array containing the collection interval, in ms
%   cax: array of CAX/Normalized dose values, in cGy
%   datatype: string containing the measurement data type
%   tics: array of timer tics, in microseconds
%   xdata: 2D array of X axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the data for each measurement
%   ydata: 2D array of Y axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the data for each measurement
%   pdiag: 2D array of positive diagonal data, where column one is the
%       position (in cm), and columns 2:n+1 are the data for each 
%       measurement
%   ndiag: 2D array of negative diagonal data, where column one is the
%       position (in cm), and columns 2:n+1 are the data for each 
%       measurement
%
% Below is an example of how this function is used:
%
%   % Load SNC ASCII data
%   path = '/path/to/files/';
%   name = 'Head1_G0.txt';
%   data = ParseSNCtxt(path, name);
%
%   % Plot Y axis profiles
%   figure;
%   hold on;
%   for i = 2:length(data.ydata)
%       plot(data.ydata{1}, data.ydata{i} * data.cax(i-1));
%   end
%   hold off;
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
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

% Define run parameter defaults
opt.Progress = true;

% Update runtime parameter arguments
for i = 1:2:nargin-2
    opt.(varargin{i}) = varargin{i+1};
end

% Log start of file load and start timer
if exist('Event', 'file') == 2
    Event(['Parsing SNC ASCII file ', name]);
    tic;
end

% Attempt to open file handle to data
fid = fopen(fullfile(path, name), 'r');

% Verify file handle is valid
if fid >= 3
    if exist('Event', 'file') == 2
        Event('Read handle successfully established');
    end
else
    if exist('Event', 'file') == 2
        Event(['Read handle not successful for ', name], 'ERROR');
    else
        error(['Read handle not successful for ', name]);
    end
end

% If a valid screen size is returned (MATLAB was run without -nodisplay)
if opt.Progress && usejava('jvm') && feature('ShowFigureWindows')
    
    % Start waitbar
    progress = waitbar(0, 'Parsing SNC ASCII file');
end

% Initialize return variable
data = struct;

% Execute in try/catch statement
try

% Retrieve the first line in the file
tline = fgetl(fid);

% Search for the Filename (suggests PROFILER)
if (isfield(opt, 'Type') && strcmpi(opt.Type, 'PROFILER')) || ...
        strcmp(tline(1:8), 'Filename')

    % Calculate number of files loaded
    n = length(regexp(tline(9:end), '\t[^\t]+'));
    if exist('Event', 'file') == 2
        Event(sprintf('ASCII file contains %i profiles', n));
    end

    % Store filenames
    data.filenames = regexp(tline(9:end), '\t([^\t]+)', 'tokens');
    
    % Declare search variables. This array specifies what lines are extracted
    % from the file, and into what format. The first column is the stored 
    % variable name, the second is the search string, and the third is the data 
    % type. Currently supported type values are string, float, datenum, and 
    % array; see while loop for specifics.
    search = {
        'timestamp'     'TimeStamp'             'datenum'
        'description'   'Description'           'string'
        'institution'   'Institution'           'string'
        'version'       'Software Version'      'string'
        'room'          'Room'                  'string'
        'mtype'         'Machine Type'          'string'
        'mmodel'        'Machine Model'         'string'
        'mserial'       'Machine Serial Number' 'string'
        'mbeamtype'     'Beam Type'             'string'
        'menergy'       'Energy'                'string'
        'wangle'        'Wedge Angle'           'float'
        'wtype'         'Wedge Type'            'string'
        'mangle'        'Gantry Angle'          'float'
        'cangle'        'Collimator Angle'      'float'
        'cleft'         'Collimator Left'       'float'
        'cright'        'Collimator Right'      'float'
        'ctop'          'Collimator Top'        'float'
        'cbottom'       'Collimator Bottom'     'float'
        'mrate'         'Rate'                  'float'
        'mdose'         'Dose'                  'float'
        'dorientation'  'Orientation'           'string'
        'dssd'          'SSD'                   'float'
        'dbuildup'      'Buildup'               'float'
        'dcal'          'Calibration File'      'string'
        'dmodel'        'Collector Model'       'string'
        'dserial'       'Collector Serial'      'string'
        'drev'          'Collector Revision'    'string'
        'dfirmware'     'Firmware Version'      'string'
        'dmode'         'Measurement Mode'      'string'
        'dgain'         'Nominal Gain'          'float'
        'dinterval'     'Collection Interval'   'float'
        'cax'           'CAX Dose'              'float'
        'datatype'      'Measured Data:'        'string'
        'tics'          'TimerTics'             'float'
        'xdata'         'Detector ID	X Axis Position(cm)'  'array'
        'ydata'         'Detector ID	Y Axis Position(cm)'  'array'
        'pdiag'         'Detector ID	Positive Diagonal Position(cm)'  'array'
        'ndiag'         'Detector ID	Negative Diagonal Position(cm)'  'array'
    };

% Otherwise, search for ** (suggests ARCCHECK)
elseif (isfield(opt, 'Type') && strcmpi(opt.Type, 'ARCCHECK')) || ...
        strcmp(tline(1:2), '**')

    % Skip empty lines until software version
    while ~feof(fid)
        tline = fgetl(fid);
        if ~isempty(tline)
            break
        end
    end
    data.version = tline;
    
    % Store file revision
    fields = strsplit(fgetl(fid), ':');
    data.filerev = strip(fields{2});
    if ~strcmp(data.filerev, 'E')
        if exist('Event', 'file') == 2
            Event('File version is incompatible with this function', ...
                'ERROR');
        else
            error('File version is incompatible with this function');
        end
    end
    
    % Store vendor, measurement device
    fgetl(fid);
    data.dvendor = strip(fgetl(fid));
    data.dmodel = strip(fgetl(fid));
    
    % Pre-define number of detectors + 1 (ArcCHECK)
    n = 132;
    
    % Declare search variables. This array specifies what lines are extracted
    % from the file, and into what format. The first column is the stored 
    % variable name, the second is the search string, and the third is the data 
    % type. Currently supported type values are string, float, datenum, and 
    % array; see while loop for specifics.
    search = {
        'filename'      'FileName:'                 'string'
        'dfirmware'     'Firmware Version:'         'string'
        'drev'          'Hardware Revision:'        'string'
        'dtype'         'Diode Type:'               'string'
        'dtemp'         'Temperature:'              'float'
        'dtilt'         'Inclinometer Tilt:'        'float'
        'drotation'     'Inclinometer Rotation:'    'float'
        'dthreshold'    'Background Threshold:'     'float'
        'dplugdose'     'Measured Cavity Dose:'     'floatunits'
        'timestamp'     'Date:'                     'dateandtime'
        'dserial'       'Serial No:'                'string'
        'doverrange'    'Overrange Error:'          'float'
        'dcal'          'Cal File:'                 'string'
        'dcalvalue'     'Dose per Count:'           'float'
        'dcalinfo'      'Dose Info:'                'string'
        'dcaliddc'      'Dose IDDC:'                'string'
        'dtime'         'Time:'                     'float'
        'dorientation'  'Orientation:'              'logical'
        'drows'         'Rows:'                     'float'
        'dcols'         'Cols:'                     'float'
        'dcaxx'         'CAX X:'                    'float'
        'dcaxy'         'CAX Y:'                    'float'
        'dqa'           'Device Position QA:'       'float'
        'shiftx'        'Shift X:'                  'floatunits'
        'shifty'        'Shift Y:'                  'floatunits'
        'shiftz'        'Shift Z:'                  'floatunits'
        'rotx'          'Rotation X:'               'floatunits'
        'rtoy'          'Rotation Y:'               'floatunits'
        'rotz'          'Rotation Z:'               'floatunits' 
        'mtype'         'Manufacturer:'             'string'
        'menergy'       'Energy:'                   'string'
        'dplug'         'Plug Present:'             'logical'
        'angular'       'Applied Angular:'          'logical'
        'fieldsize'     'Applied Field Size:'       'logical'
        'hetero'        'Applied Heterogeneity:'    'logical' 
        'background'    'Ycm	ROW'                'array'
        'calfactors'    'Ycm	ROW'                'array'
        'offset'        'Ycm	ROW'                'array'
        'rawcounts'     'Ycm	ROW'                'array'
        'corrcounts'    'Ycm	ROW'                'array'
        'dosecounts'    'Ycm	ROW'                'array'
        'flags'         'Ycm	ROW'                'array'
        'interpolated'  'Ycm	ROW'                'array'
        'doseinterp'    'Ycm	ROW'                'array'
    };
    
else
    
    % Otherwise, file may not be in correct format
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
    for i = 1:size(search, 1)
       
        % If search variable is found
        if ~isfield(data, char(search(i,1))) && ...
                length(tline) >= length(char(search(i,2)))+1 && ...
                strcmp(sprintf('%s\t', char(search(i,2))), ...
                tline(1:length(char(search(i,2)))+1))
            
            % Update waitbar
            if exist('progress', 'var') && ishandle(progress)
                waitbar(i/size(search, 1), progress);
            end
            
            % If returning a string
            if strcmp(search(i,3), 'string')
                
                % Store results as cell array of strings
                data.(char(search(i,1))) = regexp(tline(length(char(...
                    search(i,2)))+1:end), '\t{1,2}([^\t]+)', 'tokens');
                
            % Otherwise, if returning a logical
            elseif strcmp(search(i,3), 'logical')
                
                % Temporarily store cell array
                C = regexp(tline(length(char(search(i,2)))+2:end), ...
                    '\t([^\t]+)', 'tokens');
                
                % Search alternative format
                if isempty(C)
                     C = regexp(tline(length(char(search(i,2)))+2:end), ...
                        '([^,]+)', 'tokens');
                end
                
                % Loop through each result, computing logical
                for j = 1:length(C)
                    data.(char(search(i,1)))(j) = logical(str2double(C{j}));
                end
                
                % Clear temporary variables
                clear C j;
                
            % Otherwise, if returning a float
            elseif strcmp(search(i,3), 'float')
                
                % Temporarily store cell array
                C = regexp(tline(length(char(search(i,2)))+2:end), ...
                    '\t([^\t]+)', 'tokens');
                
                % Search alternative format
                if isempty(C)
                     C = regexp(tline(length(char(search(i,2)))+2:end), ...
                        '([^,]+)', 'tokens');
                end
                
                % Loop through each result, computing double
                for j = 1:length(C)
                    data.(char(search(i,1)))(j) = str2double(C{j});
                end
                
                % Clear temporary variables
                clear C j;

            % Otherwise, if returning a float with units
            elseif strcmp(search(i,3), 'floatunits')
                
                % Temporarily store cell array
                C = regexp(tline(length(char(search(i,2)))+2:end), ...
                        '([^\s]+)', 'tokens');
                
                % Store first value as double
                if ~isempty(C)
                    data.(char(search(i,1))) = str2double(C{1});
                end
                
                % Clear temporary variables
                clear C;
                
            % Otherwise, if returning a datenum
            elseif strcmp(search(i,3), 'datenum')
                
                % Temporarily store cell array
                C = regexp(tline(length(char(search(i,2)))+1:end), ...
                    '\t([^\t]+)', 'tokens');
                
                % Loop through each result, computing datenum
                for j = 1:length(C)
                    data.(char(search(i,1)))(j) = datenum(C{j});
                end
                
                % Clear temporary variables
                clear C j;
                
            % Otherwise, if returning a dateandtime
            elseif strcmp(search(i,3), 'dateandtime')
                
                % Temporarily store cell array
                C = regexp(tline(length(char(search(i,2)))+1:end), ...
                    '\t([^\t]+)', 'tokens');
                
                % Loop through each result, computing datenum
                data.(char(search(i,1))) = datenum([C{1}{1}, ' ', C{3}{1}]);
                
                % Clear temporary variables
                clear C;
                
            % Otherwise, if returning an array
            elseif strcmp(search(i,3), 'array')
                
                % Scan subsequent lines in file for array data
                data.(char(search(i,1))) = ...
                    cell2mat(textscan(fid, repmat('%f ', 1, n+1)));
                
            % Otherwise, return an error
            else
                if exist('Event', 'file') == 2
                    Event('Search variable type is not supported', 'ERROR');
                else
                    error('Search variable type is not supported');
                end
            end
            
            % End search, as a match was found
            break;
        end
    end
end

% Close file handle
fclose(fid);

% Store single cell arrays as their contents
for i = 1:size(search, 1)
    if isfield(data, search{i,1}) && iscell(data.(search{i,1})) && ...
            length(data.(search{i,1})) == 1
        data.(search{i,1}) = data.(search{i,1}){1}{1};
    end
end

% Log SNC version, if available
if exist('Event', 'file') == 2 && isfield(data, 'version') && ...
        ~isempty(data.version)
    Event(sprintf('SNC Version %s', char(data.version{1})));
end

% Log completion of function
if exist('Event', 'file') == 2
    Event(sprintf('Successfully parsed file in %0.3f seconds', toc));
end

% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

% Clean up ArcCHECK arrays
if isfield(data, 'dcols')
    for i = 1:size(search, 1)
        if isfield(data, search{i,1}) && ~ischar(data.(search{i,1})) && ...
                size(data.(search{i,1}), 2) == data.dcols + 2
            data.(search{i,1}) = data.(search{i,1})(:, 3:end);
        end
    end
end

% Clear temporary files
clear fid tline n i progress search;

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
