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
%   xfwhms: vector of X axis Full Width at Half Maximum(s) for each field
%   xedges: n x 2 array of left and right FWHM-defined X axis field edges
%   yfwhms: vector of Y axis Full Width at Half Maximum(s) for each field
%   yedges: n x 2 array of left and right FWHM-defined Y axis field edges
%   pfwhms: vector of positive diagonal axis Full Width at Half Maximum(s) 
%       for each field
%   pedges: n x 2 array of left and right FWHM-defined positive diagonal 
%       axis field edges
%   nfwhms: vector of negative diagonal axis Full Width at Half Maximum(s) 
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
% reference data was provided:
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
%   names = 'Head1_G90_27p3.prm';
%   data = ParseSNCprm(path, names);
%
%   % Compute statistics on Profiler data
%   results = AnalyzeProfilerFields(data);
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
%   [results, refresults] = AnalyzeProfilerFields(data, refdata);
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
% Initialize data field order array
fields = {'xdata', 'ydata', 'pdiag', 'ndiag'};

% Check the number of inputs
if nargin == 0 || nargin > 3
    if exist('Event', 'file') == 2
        Event('Incorrect number of arguments passed to function', 'ERROR');
    else
        error('Incorrect number of arguments passed to function');
    end
end

% Check number of input vs. output arguments
if nargout > 1 && nargin == 1
    if exist('Event', 'file') == 2
        Event(['Function cannot return two output arguments with only', ...
            ' one input argument'], 'ERROR');
    else
        error(['Function cannot return two output arguments with only', ...
            ' one input argument']);
    end
    
% Otherwise, check number of output arguments
elseif nargout == 0 || nargout > 2
    if exist('Event', 'file') == 2
        Event('Incorrect number of outputs arguments', 'ERROR');
    else
        error('Incorrect number of outputs arguments');
    end
end

% If first argument contains the field 'data'
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
elseif isstruct(varargin{1}) && isfield(varargin{1}, 'xdata')
    
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

% Execute in try/catch statement
try

% Log start of analysis and start timer
if exist('Event', 'file') == 2
    Event('Analyzing SNC Profiler data');
    tic;
end
    
%% Extract field
% If file type is PRM, split into individual fields
if strcmp(type, 'prm')
    
    % Log event
    if exist('Event', 'file') == 2
        Event('Extracting frames from time-dependent data');
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
                for j = 1:c
                    
                    % Update data array
                    varargout{1}.(fields{k})(j+1, i) = interp1(...
                        cat(2, varargout{1}.(fields{k})(1, 1:i-1), ...
                        varargout{1}.(fields{k})(1, i+1:end)), ...
                        cat(2, varargout{1}.(fields{k})(j+1, 1:i-1), ...
                        varargout{1}.(fields{k})(j+1, i+1:end)), ...
                        varargout{1}.(fields{k})(1, i), 'spline', 'extrap');
                end  
            end
        end
    end
    
    % Clear temporary variables
    clear i j k c;
    
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
    
        % Store data if present, converting from normalized dose to dose
        if isfield(varargin, fields{k})
            varargout{1}.(fields{k}) = varargin{1}.(fields{k})/100 .* ...
                repmat([1 varargin{1}.cax], ...
                size(varargin{1}.(fields{k}), 2)) / 100;
        end
    end
    
    % Clear temporary variables
    clear k;
end
    
%% Set reference profiles (if provided)  
if nargin >= 2

    % Log event
    if exist('Event', 'file') == 2
        Event('Computing reference profile correlation matrix');
    end
    
    % Initialize 3D correlation return array
    varargout{1}.corr = zeros(length(fields), size(varargout{1}.xdata, 1) ...
        - 1, size(varargin{2}.xdata, 1) - 1);
 
    % Loop through each field
    for k = 1:length(fields)

        % Initialize empty temp variable
        ref = zeros(size(varargin{2}.(fields{k}), 1), ...
            size(varargout{1}.(fields{k}), 2));
        
        % Loop through reference profiles
        for i = 2:size(varargin{2}.(fields{k}), 1)
            
            % Interpolate reference profiles to measured data using splines
            ref(i, :) = interp1(varargin{2}.(fields{k})(1, :), ...
                varargin{2}.(fields{k})(i, :), ...
                varargout{1}.(fields{k})(1, :), 'spline', 'extrap');
        end
        
        % Compute correlation coefficient of normalized profiles (note 
        % corr requires n x p1 and n x p2 arrays, so profiles must be 
        % transposed when passed)
        varargout{1}.corr(k, :, :) = corr(...
            (varargout{1}.(fields{k})(2:end, :)./...
            repmat(max(varargout{1}.(fields{k})(2:end, :), [], 2), ...
            [1, size(varargout{1}.(fields{k}), 2)]))', (ref(2:end, :)./...
            repmat(max(ref(2:end, :), [], 2), [1, size(ref, 2)]))', 'type', ...
            'Pearson');
    end
        
    % Clear temporary variables 
    clear i k ref;
    
    % Determine which reference profile has the highest correlation
    [~, varargout{1}.ref] = max(squeeze(sum(varargout{1}.corr, 1)), [], 2);    

    % Log event
    if exist('Event', 'file') == 2
        Event(sprintf(['Reference profile indices = [', repmat('%s,', ...
            [1 length(varargout{1}.ref)]), ']'], varargout{1}.ref));
    end
    
    % Set varargout{2} profile data (xdata, ydata, pdiag, ndiag)
    if nargout == 2
        
        % Loop through each field
        for k = 1:length(fields)
            
            % Store reference xdata based on highest correlation
            varargout{2}.(fields{k}) = ...
                varargin{2}.(fields{k})([1, 1 + squeeze(varargout{1}.ref)'], :);
        end
        
        % Clear temporary variables 
        clear k;
    end
end


%% Compute FWHM and field edges 
% Loop through measured and reference return variables
for i = 1:nargout
    
    % Loop through each field
    for k = 1:length(fields)

        % Find maximum value in profile
        
        
        
    end
end

% Clear temporary variables 
clear i k;
        
%% Compute flatness and symmetry
for i = 1:nargout
    



end

%% Compute profile differences and Gamma (if provided)
if nargin == 2
    

    
end
    
% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['SNC Profiler analysis successfully completed in ', ...
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


