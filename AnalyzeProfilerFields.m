function varargout = AnalyzeProfilerFields(varargin)
% AnalyzeProfilerFields reads in a data structure (created via ParseSNCprm
% or ParseSNCtxt) and computes statistics for each field.  If multiple
% frames are detected when analyzing time-dependent ParseSNCprm data
% (indicated by gaps in time between measurements), each frame will be
% separated out and statistics will be computed independently.  The
% various statistics are returned as a MATLAB structure.
%
% This function may be executed with a second input argument containing one
% or more reference profiles for each axis (xdata, ydata, pdiag, ndiag). If
% provided, additional statistics (such as Gamma comparisons) will be 
% computed. If only one reference profile is provided, but multiple
% measured profiles are provided, the same reference profile will be
% applied to each measured profile. If multiple reference profiles are
% provided, the closest matching reference profile to the measured frame 
% will be selected (using a correlation coefficient).
%
% If reference data is provided, a second structure can also be returned by
% this function containing the same fields (and array sizes) as the first
% results array, but containing the statistics for the reference profile.
%
% This function will display a progress bar while it computes Gamma (unless 
% MATLAB was executed with the -nodisplay, -nodesktop, or -noFigureWindows 
% flags).
%
% Finally, when correcting for ignored detectors or computing profile
% differences or gamma (which requires uniformly spaced data), 
% interpolation may be required.  In these instances, spline-based 
% interpolation is used.
%
% This function has been tested with SNC Profiler version 3.3.1 and .prm 
% version 25. 
%
% The following variables are required for proper execution:
%   varargin{1}: structure returned either by ParseSNCprm or ParseSNCtxt 
%       (see the documentation for these functions for more information on 
%       the fields contained)
%   varargin{2} (optional): structure containing reference profile data. 
%       The following structure fields can be incuded: xdata (2 x n+1 
%       array), ydata (2 x n+1 array), pdiag (2 x n+1 array), and ndiag 
%       (2 x n+1 array). To compute Gamma, must also include the fields
%       abs (absolute gamma criterion), dta (distance to agreement, in cm),
%       and local (flag indicating whether to perform local or global Gamma
%       comparison).
%   varargin{3} (optional): string indicating whether or not to normalize
%       profiles when processing. Can be 'none', 'max', or 'center'.  If
%       omitted, 'none' is assumed.  If reference data is not passed, this
%       may also be passed as varargin{2} (see examples).
%   varargin{4} (optional): boolean indicating whether or not to show
%       progress bar
%
% The following structure fields are returned for varargout{1} and 
% varargout{2} upon successful completion:
%   xdata: 2D array of background/ignore flag/calibration corrected X axis 
%       data, where column one is the position (in cm), and columns 2:n+1 
%       are the dose for each frame/profile (in Gy)
%   ydata: 2D array of background/ignore flag/calibration corrected Y axis 
%       data, where column one is the position (in cm), and columns 2:n+1 
%       are the dose for each frame/profile (in Gy)
%   pdiag: 2D array of background/ignore flag/calibration corrected 
%       positive diagonal axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the dose for each frame/profile (in Gy)
%   ndiag: 2D array of background/ignore flag/calibration corrected 
%       negative diagonal axis data, where column one is the position (in 
%       cm), and columns 2:n+1 are the dose for each frame/profile (in Gy)
%   tdata (optional): 2D array of time-dependent detector channel data, 
%       where column one is the absolute time and column two is the 
%       differential central channel response. If time dependent data is
%       not available (i.e. ParseSNCtxt data), this field is not returned.
%   xfwhm: vector of X axis Full Width at Half Maximum(s) for each field.
%       Note, if profile edges cannot be found for a given profile, the 
%       FWHM and edges values will be zero.
%   xedges: n x 2 array of left and right FWHM-defined X axis field edges
%   yfwhm: vector of Y axis Full Width at Half Maximum(s) for each field
%   yedges: n x 2 array of left and right FWHM-defined Y axis field edges
%   pfwhm: vector of positive diagonal axis Full Width at Half Maximum(s) 
%       for each field
%   pedges: n x 2 array of left and right FWHM-defined positive diagonal 
%       axis field edges
%   nfwhm: vector of negative diagonal axis Full Width at Half Maximum(s) 
%       for each field
%   nedges: n x 2 array of left and right FWHM-defined negative diagonal 
%       axis field edges
%   xflat: vector of central 80% flatness for X axis profiles
%   xsym: vector of areal symmetry for X axis profiles
%   yflat: vector of central 80% flatness for Y axis profiles
%   ysym: vector of areal symmetry for Y axis profiles
%   pflat: vector of central 80% flatness for positive diagonal profiles
%   psym: vector of areal symmetry for positive diagonal profiles
%   nflat: vector of central 80% flatness for negative diagonal profiles
%   nsym: vector of areal symmetry for negative diagonal profiles
%
% The following additional structure fields are returned in varargout{1} if 
% reference data was provided and nargout == 2:
%   ref (optional): vector indicating which reference profile was selected
%   corr (optional): 4 x n x m 3D array containing the correlation 
%       coefficient between each of n measured profiles, m reference 
%       profiles, and axis (x, y, p, n)
%   xdiff (optional): 2D array of X axis differences, where column one 
%       is the position (in cm), and columns 2:n+1 are the abs differences
%   ydiff (optional): 2D array of Y axis differences, where column one 
%       is the position (in cm), and columns 2:n+1 are the abs differences
%   pdiff (optional): 2D array of positive diagonal differences, where 
%       column one is the position (in cm), and columns 2:n+1 are the abs
%       differences
%   ndiff (optional): 2D array of negative diagonal differences, where 
%       column one is the position (in cm), and columns 2:n+1 are the abs
%       differences
%   xgamma (optional): 2D array of X axis Gamma values, where column one 
%       is the position (in cm), and columns 2:n+1 are the Gamma indices
%   ygamma (optional): 2D array of Y axis Gamma values, where column one 
%       is the position (in cm), and columns 2:n+1 are the Gamma indices
%   pgamma (optional): 2D array of positive diagonal Gamma values, where 
%       column one is the position (in cm), and columns 2:n+1 are the Gamma
%       indices
%   ngamma (optional): 2D array of negative diagonal Gamma values, where 
%       column one is the position (in cm), and columns 2:n+1 are the Gamma
%       indices
%
% Below is an example of how this function is used:
%
%   % Load SNC Profiler data from one file
%   path = '/path/to/files/';
%   name = 'Head1_27p3.prm';
%   data = ParseSNCprm(path, name);
%
%   % Compute statistics on Profiler data (without normalization)
%   results = AnalyzeProfilerFields(data, 'none');
%
%   % Load reference profiles from a DICOM dose file
%   refdata = LoadProfilerDICOMReference(...
%       'AP_27P3X27P3_PlaneDose_Vertical_Isocenter.dcm');
%
%   % Set Gamma criteria
%   refdata.abs = 2;    % percent
%   refdata.dta = 0.1;  % cm
%   refdata.local = 0;  % perform global analysis
%
%   % Compute statistics again, this time comparing to reference criteria
%   % and normalizing profiles to Profiler center
%   [results, refresults] = AnalyzeProfilerFields(data, refdata, 'center');
%
%   % Plot X axis Gamma index
%   figure;
%   subplot(2,2,1);
%   hold on;
%   plot(results.xdata(1,:), results.xdata(2,:));
%   plot(refresults.xdata(1,:), refresults.xdata(2,:));
%   plot(results.xgamma(1,:), results.xgamma(2,:));
%   hold off;
%   title(sprintf('xgamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
%   xlabel('x axis (cm)');
% 
%   % Plot Y axis Gamma index
%   subplot(2,2,2);
%   hold on;
%   plot(results.ydata(1,:), results.ydata(2,:));
%   plot(refresults.ydata(1,:), refresults.ydata(2,:));
%   plot(results.ygamma(1,:), results.ygamma(2,:));
%   hold off;
%   title(sprintf('ygamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
%   xlabel('y axis (cm)');
% 
%   % Plot positive diagonal axis Gamma index
%   subplot(2,2,3);
%   hold on;
%   plot(results.pdiag(1,:), results.pdiag(2,:));
%   plot(refresults.pdiag(1,:), refresults.pdiag(2,:));
%   plot(results.pgamma(1,:), results.pgamma(2,:));
%   hold off;
%   title(sprintf('pgamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
%   xlabel('positive diagonal (cm)');
% 
%   % Plot negative diagonal axis Gamma index
%   subplot(2,2,4);
%   hold on;
%   plot(results.ndiag(1,:), results.ndiag(2,:));
%   plot(refresults.ndiag(1,:), refresults.ndiag(2,:));
%   plot(results.ngamma(1,:), results.ngamma(2,:));
%   hold off;
%   title(sprintf('ngamma (%0.1f%%/%0.1f cm)', refdata.abs, refdata.dta));
%   xlabel('negative diagonal (cm)');
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2015-2018 University of Wisconsin Board of Regents
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
% Initialize data field order array
fields = {'xdata', 'ydata', 'pdiag', 'ndiag'};

% If too few or too many input arguments are provided
if nargin == 0 || nargin > 3
    
    % Throw an error
    if exist('Event', 'file') == 2
        Event('Incorrect number of arguments passed to function', 'ERROR');
    else
        error('Incorrect number of arguments passed to function');
    end
end

% If an incorrect number of input vs. output arguments exist
if (nargout > 1 && nargin == 1) || (nargout > 1 && nargin == 2 && ...
        ~isstruct(varargin{2})) 
    
    % Throw an error
    if exist('Event', 'file') == 2
        Event(['Function cannot return two output arguments with only', ...
            ' one data input argument'], 'ERROR');
    else
        error(['Function cannot return two output arguments with only', ...
            ' one data input argument']);
    end
    
% Otherwise, if too many or too few output arguments exist
elseif nargout == 0 || nargout > 2
    
    % Throw an error
    if exist('Event', 'file') == 2
        Event('Incorrect number of outputs arguments', 'ERROR');
    else
        error('Incorrect number of outputs arguments');
    end
end

% If first argument contains the correct fields
if isstruct(varargin{1}) && isfield(varargin{1}, 'data') && ...
        isfield(varargin{1}, 'num') && ...
        isfield(varargin{1}, 'detspacing') && ...
        isfield(varargin{1}, 'background') && ...
        isfield(varargin{1}, 'calibration') && ...
        isfield(varargin{1}, 'ignore')
        
    % Set type to PRM (from ParseSNCprm)
    type = 'prm';
    
    % Log event
    if exist('Event', 'file') == 2
        Event('PRM profiler data detected for analysis');
    end
    
    % Flip X & Y num values (they appear to be flipped in PRM files
    if exist('Event', 'file') == 2
        Event('Swapping X & Y number of elements');
    end
    varargin{1}.num(1:2) = varargin{1}.num(2:-1:1);
    
% Otherwise, if the first argument contains the field 'xdata'
elseif isstruct(varargin{1}) && isfield(varargin{1}, 'xdata') && ...
        isfield(varargin{1}, 'datatype') && isfield(varargin{1}, 'cax')
    
    % Set type to ASCII (from ParseSNCtxt)
    type = 'ascii';
    
    % Log event
    if exist('Event', 'file') == 2
        Event('ASCII profiler data detected for analysis');
    end

% Otherwise, type is invalid so throw an error
else
    if exist('Event', 'file') == 2
        Event('Invalid first input argument', 'ERROR');
    else
        error('Invalid first input argument');
    end
end

% Check for normalization input argument
if nargin >= 2 && ischar(varargin{2})
    
    % Set from second input
    norm = varargin{2};

% Otherwise, if three arguments are provided
elseif nargin == 3 && ischar(varargin{3})
    
    % Set from second input
    norm = varargin{3};
    
% Otherwise, no normalization value was provided
else
    
    % Default to none
    norm = 'none'; 
end

% Determine path of current application
[path, ~, ~] = fileparts(mfilename('fullpath'));

% Add gamma submodule to search path
addpath(fullfile(path, '/gamma'));

% Clear temporary variable
clear path;

% Check if MATLAB can find CalcGamma
if exist('CalcGamma', 'file') ~= 2
    
    % If not, throw an error
    if exist('Event', 'file') == 2
        Event(['The CalcGamma submodule does not exist in the search path.', ...
            ' Use  git clone --recursive or git submodule init followed by', ...
            ' git submodule update to fetch all submodules'], 'ERROR');
    else
        error(['The CalcGamma submodule does not exist in the search path.', ...
            ' Use  git clone --recursive or git submodule init followed by', ...
            ' git submodule update to fetch all submodules']);
    end
end

% Check if MATLAB can find corr (Statistics Toolbox)
if exist('corr', 'file') ~= 6 && exist('corr', 'file') ~= 2
    
    % If not, throw an error
    if exist('Event', 'file') == 2
        Event(['The Statistics and Machine Learning Toolbox cannot be ', ...
            'found and is required by this function.'], 'ERROR');
    else
        error(['The Statistics and Machine Learning Toolbox cannot be ', ...
            'found and is required by this function.']);
    end
end

% Execute in try/catch statement
try

% Log start of analysis and start timer
if exist('Event', 'file') == 2
    Event('Analyzing SNC Profiler data');
    timer = tic;
end 

% Verify normalization value
if strcmp(norm, 'none') || strcmp(norm, 'center') || strcmp(norm, 'max')
   
    % Log value
    if exist('Event', 'file') == 2
        Event(['Normalization value set to ', norm]);
    end
    
% Otherwise, throw an error    
else    
    if exist('Event', 'file') == 2
        Event('Normalization input parameter is invalid', 'ERROR');
    else
        error('Normalization input parameter is invalid');
    end
end

%% Extract field
% If file type is PRM, split into individual fields
if strcmp(type, 'prm')
    
    % If the dosecal value is 0, warn user and override to 1
    if varargin{1}.dosecal == 0
        if exist('Event', 'file') == 2
            Event('Dose calibration factor is not set, ignoring', 'WARN');
        end
        varargin{1}.dosecal = 1;
    end
    
    % Log event
    if exist('Event', 'file') == 2
        Event('Initializing detector spatial positions');
    end
      
    % Initialize X axis positions (note that two diodes are missing)
    varargout{1}.xdata(1, 1:varargin{1}.num(1)) = ...
        ([1:(varargin{1}.num(1) - 1) / 2, (varargin{1}.num(1) + 3) / 2, ...
        (varargin{1}.num(1) + 7) / 2:varargin{1}.num(1) + 2] - ...
        (varargin{1}.num(1) + 3) / 2) * varargin{1}.detspacing;

    % Initialize Y axis positions
    varargout{1}.ydata(1, 1:varargin{1}.num(2)) = ((1:varargin{1}.num(2)) - ...
        (varargin{1}.num(2) + 1) / 2) * varargin{1}.detspacing;
    
    % Initialize pos diagonal position (note that two diodes are missing)
    varargout{1}.pdiag(1, 1:varargin{1}.num(3)) = ...
        ([1:(varargin{1}.num(3) - 1) / 2, (varargin{1}.num(3) + 3) / 2, ...
        (varargin{1}.num(3) + 7) / 2:varargin{1}.num(3) + 2] - ...
        (varargin{1}.num(3) + 3) / 2) * varargin{1}.detspacing * sqrt(2);

    % Initialize neg diagonal position (note that two diodes are missing)
    varargout{1}.ndiag(1, 1:varargin{1}.num(4)) = ...
        ([1:(varargin{1}.num(4) - 1) / 2, (varargin{1}.num(4) + 3) / 2, ...
        (varargin{1}.num(4) + 7) / 2:varargin{1}.num(4) + 2] - ...
        (varargin{1}.num(4) + 3) / 2) * varargin{1}.detspacing * sqrt(2);
    
    % Log event
    if exist('Event', 'file') == 2
        Event('Extracting frames from time-dependent data');
    end
    
    % Initialize counters
    c = 0;
    i = 0;
    
    % Continue while the counter has not reached the end of the data array
    while i < size(varargin{1}.data, 1)
        
        % Increment data point counter
        i = i + 1;
    
        % Begin searching ahead
        for j = i+1:size(varargin{1}.data, 1)

            % If there is a gap in the data
            if j == size(varargin{1}.data, 1) || ...
                    abs(varargin{1}.data(j+1, 2) - ...
                    varargin{1}.data(j, 2)) > 2 

                % Increment frame counter
                c = c + 1;
                
                % Log find
                if exist('Event', 'file') == 2
                    Event(sprintf('Frame %i identified from %i to %i', ...
                        c, i, j));
                end
                
                % Loop through the data fields
                for k = 1:min(length(fields), length(varargin{1}.num))
                
                    % Store field data, correcting for background/cal/dose,
                    % skipping the first 5 data columns
                    varargout{1}.(fields{k})(c+1, 1:varargin{1}.num(k)) = ...
                        ((varargin{1}.data(j, 5 + ...
                        sum(varargin{1}.num(1:k-1)) + ...
                        (1:varargin{1}.num(k))) - varargin{1}.data(i, 5 + ...
                        sum(varargin{1}.num(1:k-1)) + ...
                        (1:varargin{1}.num(k)))) - ((varargin{1}.data(j, 3) ...
                        - varargin{1}.data(i, 3)) * varargin{1}.background(...
                        (2 + sum(varargin{1}.num(1:k-1))):(1 + ...
                        varargin{1}.num(k) + ...
                        sum(varargin{1}.num(1:k-1)))))) .* ...
                        varargin{1}.calibration((1 + ...
                        sum(varargin{1}.num(1:k-1))):(varargin{1}.num(k) + ...
                        sum(varargin{1}.num(1:k-1)))) * ...
                        varargin{1}.dosecal / 1000;
                end
                
                % Jump forward
                i = j;
                break;      
            end
        end
    end
    
    % Log event
    if exist('Event', 'file') == 2
        Event(sprintf('Extraction complete, %i frames loaded', c));
    end
    
    % Clear temporary variables
    clear i j k;

    % Log event
    if exist('Event', 'file') == 2
        Event('Checking ignore flags');
    end
    
    % Loop through the data fields
    for k = 1:min(length(fields), length(varargin{1}.num))
        
        % Loop through the detectors in field k
        for i = 1:varargin{1}.num(k)
            
            % If ignore detector flag is set
            if varargin{1}.ignore(sum(varargin{1}.num(1:k-1)) + i) == 1
                
                % Log event
                if exist('Event', 'file') == 2
                    Event(sprintf(['%s detector %i ignored, interpolated', ...
                    ' from neighboring values'], fields{k}, i));
                end
                
                % Interpolate each ignored detector using splines
                varargout{1}.(fields{k})(2:end, i) = interp1(...
                    cat(2, varargout{1}.(fields{k})(1, 1:i-1), ...
                    varargout{1}.(fields{k})(1, i+1:end))', ...
                    cat(2, varargout{1}.(fields{k})(2:end, 1:i-1), ...
                    varargout{1}.(fields{k})(2:end, i+1:end))', ...
                    varargout{1}.(fields{k})(1, i)', 'spline', 'extrap')';
            end
        end
    end
    
    % Clear temporary variables
    clear i k c;
  
    % Load time-dependent reference detector profile
    if exist('Event', 'file') == 2
        Event('Loading center detector timing profile');
    end
    
    % Store time values in column 1 using collection interval (in msec)
    varargout{1}.tdata(1,:) = (1:size(varargin{1}.data, 1)) * ...
        varargin{1}.dinterval;
    
    % Store cumulative raw central detector response (central Y detector)
    varargout{1}.tdata(2,:) = varargin{1}.data(:, 5 + ...
        varargin{1}.num(1) + ceil(varargin{1}.num(2)/2));
    
    % Convert signal from integral to differential using circshift
    if exist('Event', 'file') == 2
        Event('Converting cumulative timing profile to differential');
    end
    varargout{1}.tdata(2,:) = varargout{1}.tdata(2,:) - ...
        circshift(varargout{1}.tdata(2,:),1,2);
    
    % Fix first value (artifact of using circshift)
    varargout{1}.tdata(2,1) = varargout{1}.tdata(2,2);
    
% If file type is ASCII, just copy from input data
elseif strcmp(type, 'ascii')
    
    % Log event
    if exist('Event', 'file') == 2
        Event('Loading ASCII profiles');
    end
    
    % Loop through the data fields
    for k = 1:length(fields)
    
        % If the data exists
        if isfield(varargin{1}, fields{k})
            
            % Increase num variable
            varargin{1}.num(k) = size(varargin{1}.(fields{k}), 1);
            
            % If data type is normalized
            if strcmp(varargin{1}.datatype{1}, 'Total Dose-Normalized (%)')
                
                % Store data, converting from normalized dose to dose (Gy)
                varargout{1}.(fields{k}) = varargin{1}.(fields{k})'/100 .* ...
                    repmat([100;varargin{1}.cax'/100], ...
                    [1 size(varargin{1}.(fields{k}), 1)]);
            
            % Otherwise, if data type is total dose
            elseif isempty(regexp(varargin{1}.datatype{1}, 'Total Dose', ...
                    'Once'))
                
                % Store data, converting from cGy to Gy
                varargout{1}.(fields{k}) = varargin{1}.(fields{k})/100;
            
            % Otherwise, throw an error
            else
                if exist('Event', 'file') == 2
                    Event('ASCII data type is unknown', 'ERROR');
                else
                    error('ASCII data type is unknown');
                end
            end
        else
            
            % Otherwise stop loading, as not all profile data was found
            break;
        end
    end
    
    % Clear temporary variables
    clear k;
end
    
%% Set reference profiles (if provided)  
if nargin >= 2 && isstruct(varargin{2})

    % Log event
    if exist('Event', 'file') == 2
        Event('Computing reference profile correlation matrix');
    end
    
    % Initialize 3D correlation return array
    varargout{1}.corr = zeros(length(fields), size(varargout{1}.xdata, 1) ...
        - 1, size(varargin{2}.xdata, 1) - 1);
 
    % Loop through each field
    for k = 1:min(length(fields), length(varargin{1}.num))

        % Interpolate reference profiles to measured data using splines
        ref = interp1(varargin{2}.(fields{k})(1, :)', ...
            varargin{2}.(fields{k})(2:end, :)', ...
            varargout{1}.(fields{k})(1, :)', 'spline', 'extrap')';

        % Compute correlation coefficient of normalized profiles (note 
        % corr requires n x p1 and n x p2 arrays, so profiles must be 
        % transposed when passed)
        varargout{1}.corr(k, :, :) = corr(...
            (varargout{1}.(fields{k})(2:end, :)./...
            repmat(max(varargout{1}.(fields{k})(2:end, :), [], 2), ...
            [1, size(varargout{1}.(fields{k}), 2)]))', (ref ./ ...
            repmat(max(ref, [], 2), [1, size(ref, 2)]))', 'type', ...
            'Pearson');
    end
        
    % Clear temporary variables 
    clear k ref;
    
    % Determine which reference profile has the highest correlation
    if size(varargout{1}.corr, 2) == 1
        [~, varargout{1}.ref] = ...
            max(squeeze(sum(varargout{1}.corr, 1)), [], 1);
    else
        [~, varargout{1}.ref] = ...
            max(squeeze(sum(varargout{1}.corr, 1)), [], 2);
    end

    % Log event
    if exist('Event', 'file') == 2
        Event(['Reference profile indices = [', strjoin(cellstr(...
            int2str(varargout{1}.ref)), ', '), ']']);
    end
    
    % Set varargout{2} profile data (xdata, ydata, pdiag, ndiag)
    if nargout == 2
        
        % Loop through each field
        for k = 1:min(length(fields), length(varargin{1}.num))
            
            % Store reference xdata based on highest correlation
            varargout{2}.(fields{k}) = varargin{2}.(fields{k})([1, ...
                1 + squeeze(varargout{1}.ref)'], :);
        end
        
        % Clear temporary variables 
        clear k;
    end
end

%% Normalize data (optional)
% If normalization is set to center
if strcmp(norm, 'center')

    % Log event
    if exist('Event', 'file') == 2
        Event('Normalizing profiles to center value');
    end

    % Loop through measured and reference return variables
    for i = 1:nargout
        
        % Loop through the data fields
        for k = 1:min(length(fields), length(varargin{1}.num))
    
            % Normalize fields by center value
            varargout{i}.(fields{k})(2:end, :) = ...
                varargout{i}.(fields{k})(2:end, :) ./ ...
                repmat(interp1(varargout{i}.(fields{k})(1,:)', ...
                varargout{i}.(fields{k})(2:end,:)', ...
                0, 'spline')', [1 size(varargout{i}.(fields{k}),2)]);
        end
    end
    
% Otherwise, if set to max
elseif strcmp(norm, 'max')

    % Log event
    if exist('Event', 'file') == 2
        Event('Normalizing profiles to maximum value');
    end
    
    % Loop through measured and reference return variables
    for i = 1:nargout
        
        % Loop through the data fields
        for k = 1:min(length(fields), length(varargin{1}.num))
    
            % Normalize fields by max value
            varargout{i}.(fields{k})(2:end, :) = ...
                varargout{i}.(fields{k})(2:end, :) ./ ...
                repmat(max(varargout{i}.(fields{k})(2:end,:), [], 2), ...
                [1 size(varargout{i}.(fields{k}),2)]);
        end
    end
end

%% Compute FWHM and field edges 
if exist('Event', 'file') == 2
    Event('Computing full width at half maximum');
end

% Loop through measured and reference return variables
for i = 1:nargout
    
    % Loop through each field
    for k = 1:min(length(fields), length(varargin{1}.num))

        % Initialize edges return structure fields
        varargout{i}.([fields{k}(1), 'edges']) = ...
            zeros(size(varargout{i}.(fields{k}), 1) - 1, 2); %#ok<*AGROW>
        
        % Find maximum value in profile
        [C, I] = max(varargout{i}.(fields{k}), [], 2);
        
        % Loop through each profile
        for j = 2:size(varargout{i}.(fields{k}), 1)
           
            % Find highest lower index just below half maximum
            lI = find(varargout{i}.(fields{k})(j, ...
                1:I(j)) < C(j)/2, 1, 'last');

            % Find lowest upper index just above half maximum
            uI = find(varargout{i}.(fields{k})(j, ...
                I(j):end) < C(j)/2, 1, 'first');
            
            % Verify edges were found
            if isempty(uI) || isempty(lI)
                
                % Log event
                if exist('Event', 'file') == 2
                    Event(sprintf(['Field edges were not found for ', ...
                        'profile %i %s'], j-1, fields{k}), 'WARN');
                end
                
            % Otherwise, verify edges are sufficiently far from array edges
            elseif lI-1 < 1 || lI+2 > size(varargout{i}.(fields{k}), 2) || ...
                    I(j)+uI-3 < 1 || I(j)+uI > ...
                    size(varargout{i}.(fields{k}), 2)
                
                % Log event
                if exist('Event', 'file') == 2
                    Event(sprintf(['Field edges are too close to detector', ...
                        ' edge to compute for profile %i %s'], j-1, ...
                        fields{k}), 'WARN');
                end
            
            % Otherwise, continue edge calculation
            else
                
                % Interpolate to find lower half-maximum value
                varargout{i}.([fields{k}(1), 'edges'])(j-1, 1) = ...
                    interp1(varargout{i}.(fields{k})(j, lI-1:lI+2), ...
                    varargout{i}.(fields{k})(1, lI-1:lI+2), C(j)/2, ...
                    'spline');

                % Interpolate to find upper half-maximum value
                varargout{i}.([fields{k}(1), 'edges'])(j-1, 2) = ...
                    interp1(varargout{i}.(fields{k})(j, I(j)+uI-3:I(j)+uI), ...
                    varargout{i}.(fields{k})(1, I(j)+uI-3:I(j)+uI), C(j)/2, ...
                    'spline');
            end
            
            % Clear temporary variables
            clear lI uI;
        end
        
        % Compute FWHM
        varargout{i}.([fields{k}(1), 'fwhm']) = ...
            abs(varargout{i}.([fields{k}(1), 'edges'])(:,2) - ...
            varargout{i}.([fields{k}(1), 'edges'])(:,1));
        
        % Clear temporary variables
        clear C I;
    end
end

% Clear temporary variables 
clear i j k;
        
%% Compute flatness and symmetry
% Declare threshold for computing flatness and symmetry
t = 0.8;

% Set area interpolation factor (for computing area under field)
a = 10000;

% Log event
if exist('Event', 'file') == 2
    Event(sprintf('Computing central %i%% flatness and areal symmetry', ...
        t * 100));
end

% Loop through measured and reference return variables
for i = 1:nargout
    
    % Loop through each field
    for k = 1:min(length(fields), length(varargin{1}.num))
        
        % Initialize flatness and symmetry return structure fields
        varargout{i}.([fields{k}(1), 'flat']) = ...
            zeros(size(varargout{i}.(fields{k}), 1) - 1, 1);
        varargout{i}.([fields{k}(1), 'sym']) = ...
            zeros(size(varargout{i}.(fields{k}), 1) - 1, 1);
    
        % Find maximum value in profile
        [C, I] = max(varargout{i}.(fields{k}), [], 2);
    
        % Loop through each profile
        for j = 2:size(varargout{i}.(fields{k}), 1)
           
            % Find highest lower index just below half maximum
            lI = find(varargout{i}.(fields{k})(j, ...
                1:I(j)) < C(j)/2, 1, 'last');
            
            % Find lowest upper index just above half maximum
            uI = find(varargout{i}.(fields{k})(j, ...
                I(j):end) < C(j)/2, 1, 'first');
            
            % Verify edges were found
            if isempty(uI) || isempty(lI)
            
                % Log event
                if exist('Event', 'file') == 2
                    Event(sprintf(['Flatness and symmetry not computed for', ...
                        ' profile %i %s'], j-1, fields{k}), 'WARN');
                end
            
            % Otherwise, continue flatness/symmetric calculations
            else
                
                % Find maximum and minimum values within threshold region
                dmax = max(varargout{i}.(fields{k})(j, ceil(lI + (I(j) + ...
                    uI - lI) * (1 - t) / 2):floor(I(j) - 2 + uI - (I(j) + ...
                    uI - lI) * (1 - t) / 2)));
                dmin = min(varargout{i}.(fields{k})(j, ceil(lI + (I(j) + ...
                    uI - lI) * (1 - t) / 2):floor(I(j) - 2 + uI - (I(j) + ...
                    uI - lI) * (1 - t) / 2)));

                % Store flatness
                varargout{i}.([fields{k}(1), 'flat'])(j-1) = ...
                    (dmax - dmin) / (dmax + dmin);

                % Compute lower area
                lA = interp1(varargout{i}.(fields{k})(1,:), ...
                    varargout{i}.(fields{k})(j,:), varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 1) + (varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 2) - varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 1)) * (1 - t)/2:((...
                    varargout{i}.([fields{k}(1), 'edges'])(j-1, 1) + ...
                    varargout{i}.([fields{k}(1), 'edges'])(j-1, 2))/2 - ...
                    (varargout{i}.([fields{k}(1), 'edges'])(j-1, 1) + ...
                    (varargout{i}.([fields{k}(1), 'edges'])(j-1, 2) - ...
                    varargout{i}.([fields{k}(1), 'edges'])(j-1, 1)) * ...
                    (1 - t)/2)) / a:(varargout{i}.([fields{k}(1), ...
                    'edges'])(j-1, 1) + varargout{i}.([fields{k}(1), ...
                    'edges'])(j-1, 2))/2);

                % Compute upper area
                uA = interp1(varargout{i}.(fields{k})(1,:), ...
                    varargout{i}.(fields{k})(j,:), (varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 1) + varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 2))/2:((varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 2) - (varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 2) - varargout{i}.(...
                    [fields{k}(1), 'edges'])(j-1, 1)) * (1 - t)/2) - ...
                    (varargout{i}.([fields{k}(1), 'edges'])(j-1, 1) + ...
                    varargout{i}.([fields{k}(1), 'edges'])(j-1, 2))/2) / a:...
                    varargout{i}.([fields{k}(1), 'edges'])(j-1, 2) - ...
                    (varargout{i}.([fields{k}(1), 'edges'])(j-1, 2) - ...
                    varargout{i}.([fields{k}(1), 'edges'])(j-1, 1)) * ...
                    (1 - t)/2);

                % Store symmetry
                varargout{i}.([fields{k}(1), 'sym'])(j-1) = ...
                    (sum(uA) - sum(lA)) / (sum(uA) + sum(lA)) * 2;
            end
            
            % Clear temporary variables
            clear lI uI dmax dmin lA uA;
        end
        
        % Clear temporary variables
        clear C I;
    end
end

% Clear temporary variables 
clear i j k t a;

%% Compute profile differences and Gamma (if provided)
if nargout == 2 && isfield(varargin{2}, 'abs') && isfield(varargin{2}, 'dta')

    % If a valid screen size is returned (MATLAB was run without -nodisplay)
    if usejava('jvm') && feature('ShowFigureWindows') && ...
            (nargin < 4 || varargin{4})

        % Start waitbar
        progress = waitbar(0, 'Analyzing SNC Profiler data');
    end
    
    % Log event
    if exist('Event', 'file') == 2
        Event(sprintf(['Computing profile differences and %0.1f%%/%0.1fcm', ...
            ' Gamma index'], varargin{2}.abs, varargin{2}.dta));
    end
    
    % If a local gamma variable does not exist, assume global
    if ~isfield(varargin{2}, 'local')
        varargin{2}.local = 0;
        
        % Log event
        if exist('Event', 'file') == 2
            Event('Gamma index assumed to be global');
        end
    end
    
    % Loop through each field
    for k = 1:min(length(fields), length(varargin{1}.num))
    
        % Set difference and Gamma profile positions
        varargout{1}.([fields{k}(1), 'diff'])(1,:) = ...
            varargout{1}.(fields{k})(1,:);      
        varargout{1}.([fields{k}(1), 'gamma'])(1,:) = ...
            varargout{1}.(fields{k})(1,:);
        
        % Resample reference profile to measured detector positions and
        % compute difference using splines
        varargout{1}.([fields{k}(1), 'diff'])(2:size(...
            varargout{1}.(fields{k}), 1),:) = ...
            varargout{1}.(fields{k})(2:size(...
            varargout{1}.(fields{k}), 1),:) - ...
            interp1(varargout{2}.(fields{k})(1,:)', ...
            varargout{2}.(fields{k})(2:end,:)', ...
            varargout{1}.(fields{k})(1,:)', 'spline', 'extrap')';
        
        % Loop through each profile
        for j = 2:size(varargout{1}.(fields{k}), 1)

            % Log event
            if exist('Event', 'file') == 2
                Event(sprintf('Computing profile %i %sgamma', j-1, ...
                    fields{k}(1)));
            end
            
            % Update waitbar
            if exist('progress', 'var') && ishandle(progress)
                waitbar(((size(varargout{1}.(fields{k}), 1)-1) *...
                    (k-1) + (j-1))/(min(length(fields), ...
                    length(varargin{1}.num)) * ...
                    (size(varargout{1}.(fields{k}), 1)-1)), progress);
            end
            
            % Prepare CalcGamma inputs (which uses start/width/data format)
            tar.start = varargout{1}.(fields{k})(1,1);
            tar.width = varargout{1}.(fields{k})(1,2) - tar.start;
            tar.data = interp1(varargout{1}.(fields{k})(1,:), ...
                varargout{1}.(fields{k})(j,:), tar.start:tar.width:...
                varargout{1}.(fields{k})(1,end), 'spline', 'extrap');

            ref.start = varargout{2}.(fields{k})(1,1);
            ref.width = varargout{2}.(fields{k})(1,2) - ref.start;
            ref.data = interp1(varargout{2}.(fields{k})(1,:), ...
                varargout{2}.(fields{k})(j,:), ref.start:ref.width:...
                varargout{2}.(fields{k})(1,end), 'spline', 'extrap');
            
            % Calculate 1-D Gamma using CPU
            gamma = CalcGamma(ref, tar, varargin{2}.abs, varargin{2}.dta, ...
                'local', varargin{2}.local, 'cpu', 1);
            
            % Retrieve Gamma values for each detector position
            varargout{1}.([fields{k}(1), 'gamma'])(j,:) = interp1(...
                tar.start:tar.width:varargout{1}.(fields{k})(1,end), gamma, ...
                varargout{1}.(fields{k})(1,:), 'nearest', 0);
            
            % Clear temporary variables
            clear tar ref;
        end
    end
    
    % Close waitbar
    if exist('progress', 'var') && ishandle(progress)
        close(progress);
    end

    % Clear temporary variables
    clear progress;
end

% Clear temporary variables 
clear j k;
    
% Log completion
if exist('Event', 'file') == 2
    
    % Log event 
    Event(sprintf(['SNC Profiler analysis successfully completed in ', ...
      '%0.3f seconds'], toc(timer)));
  
    % Clear temporary variable
    clear timer;
end

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
