function results = AnalyzeACFields(varargin)
% AnalyzeACFields reads in an ArcCHECK data structure (created via 
% ParseSNCacm), identifies individual static fields, and computes the 
% FWHM-defined center of each field (in cylindrical coordinates). Each 
% field detector data is interpolated to a 2D cylindrical frame. The 
% results are returned as a MATLAB structure.
%
% This function will display a progress bar while it loads (unless MATLAB 
% was executed with the -nodisplay, -nodesktop, or -noFigureWindows flags).
%
% This structure uses a cylindrical coordinate system to define ArcCHECK 
% data, where Y is positioned along the central long axis of the ArcCHECK
% and theta is the angle (in degrees) from 0, where 0 is the ArcCHECK zero
% position (typically positioned in the positive IEC-Z direction)
%
% This function has been tested with SNC Patient version 6.2.3 and .acm 
% file rev C. 
% 
% When computing the 2D cylindrical frame, the MATLAB scatteredInterpolant
% class is used with a bilinear method.  When computing the FWHM-defined
% edges, splines-based interpolation is applied.
%
% The following variables are required for proper execution:
%   varargin{1}: structure returned either by ParseSNCacm (see the 
%       documentation for this function for more information on the fields 
%       contained)
%   varargin{2} (optional): number indicating the angle (in degrees) to 
%       offset the measured angles, to account for ArcCheck roll
%
% The following structure fields are returned upon successful completion:
%   detectors: 1386 x n array of detector data, corrected for background,
%       relative calibration, and absolute calibration (in Gy)
%   itheta: 360 x 201 meshgrid of theta values
%   iY: 360 x 201 meshgrid of y values
%   frames: 360 x 201 x n array of n fields, where each 2D array is an
%   	interpolated cylindrical map (theta x Y) of measured dose (in Gy)
%   alpha: 2 x n array of entrance and exit FWHM-defined profile center
%       theta values (in degrees). Note, if profile edges cannot be found 
%       for a given profile, the alpha and beta values will be zero.
%   beta: 2 x n array of entrance and exit FWHM-defined profile center Y 
%       axis values (in cm)
%
% Below is an example of how this function is used:
%
%   % Load SNC ArcCHECK data from one file
%   path = '/path/to/files/';
%   names = 'Head3_G270_to_G90_10deg.acm';
%   data = ParseSNCacm(path, names);
%
%   % Analyze ArcCHECK data
%   results = AnalyzeACFields(data);
%
%   % Loop through frames, plotting result
%   figure;
%   for i = 1:size(results.frames, 3)
%
%       % Plot frame
%       imagesc(circshift(results.frames(:,:,i), -180, 2));
%       set(gca,'XTick', 1:30:361);
%       set(gca,'XTickLabel', -180:30:180);
%       xlabel('ArcCHECK Angle (deg)');
%       set(gca,'YTick', 1:20:201);
%       set(gca,'YTickLabel', 10:-2:-10);
%       ylabel('ArcCHECK Y (cm)');
%       title(sprintf('Frame %i', i));
%   
%       % Update plot and pause temporarily
%       drawnow;
%       pause(0.1);
%   end
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

%% Validate inputs
% If too few or too many input arguments are provided
if nargin == 0 || nargin > 2
    
    % Throw an error
    if exist('Event', 'file') == 2
        Event('Incorrect number of arguments passed to function', 'ERROR');
    else
        error('Incorrect number of arguments passed to function');
    end
end

% If first argument does not contain the correct fields
if ~(isstruct(varargin{1}) && isfield(varargin{1}, 'data') && ...
        isfield(varargin{1}, 'y') && ...
        isfield(varargin{1}, 'theta') && ...
        isfield(varargin{1}, 'background') && ...
        isfield(varargin{1}, 'calibration'))
        
    % Throw an error
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

% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    
    % Start waitbar
    progress = waitbar(0, 'Analyzing SNC ArcCHECK data');
end

% Calculate interpolation meshgrids
[results.itheta, results.iY] = meshgrid(0:359, -10:0.1:10);
[~, center] = min(abs(results.iY(:,1)));

% Log interpolation grid settings
if exist('Event', 'file') == 2
    Event('Initializing detector spatial grid [0:359 deg, -10:0.1:10 cm]');
end

%% Extract frames
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
        if j == size(varargin{1}.data, 1) || abs(varargin{1}.data(j+1, 2) - ...
                varargin{1}.data(j, 2)) > 2 

            % Increment frame counter
            c = c + 1;

            % Log find
            if exist('Event', 'file') == 2
                Event(sprintf('Frame %i identified from %i to %i', ...
                    c, i, j));
            end

            % Extract data
            results.detectors(1:1386, c) = (varargin{1}.data(j, 12:1397) - ...
                varargin{1}.data(i, 12:1397) - ...
                (varargin{1}.data(j,3) - ...
                varargin{1}.data(i,3)) .* ...
                varargin{1}.background(2:1387)) .* ...
                varargin{1}.calibration(2:1387) * varargin{1}.dosecal/100;

            % Jump forward
            i = j;
            break;      
        end
    end
end

% Clear temporary variables 
clear i j;

%% Interpolate detectors into 2D frames
% Initialize frames return array
results.frames = zeros(size(results.itheta, 1), ...
    size(results.itheta, 2), c);

% Log event
if exist('Event', 'file') == 2
    Event('Computing 2D cylindrical frames');
end

% Loop through frames
for i = 1:c
    
    % Update waitbar
    if exist('progress', 'var') && ishandle(progress)
        waitbar(i/(c+1), progress);
    end

    % Log action
    if exist('Event', 'file') == 2
        Event(sprintf(['Interpolating diode data to cylindrical ', ...
            'map for frame %i'], i));
    end
    
    % Generate scattered interpolant object for diode data
    scatter = scatteredInterpolant(double([varargin{1}.theta-360, ...
        varargin{1}.theta, varargin{1}.theta+360]'), ...
        double([varargin{1}.y, varargin{1}.y, varargin{1}.y]'), ...
        double([results.detectors(:, i)', results.detectors(:, i)', ...
        results.detectors(:, i)']'), 'linear', 'linear');

    % Interpolate group data into 2D array of itheta, iY
    results.frames(:, :, i) = scatter(results.itheta, results.iY);
    
    % Clear temporary variables 
    clear scatter;
end

% Clear temporary variables
clear i;

%% Compute alpha and beta values
% Initialize alpha and beta return arrays
results.alpha = zeros(2, c);
results.beta = zeros(2, c);

% Log event
if exist('Event', 'file') == 2
    Event('Computing frame FWHM-defined center coordinates');
end

% Loop through frames
for i = 1:c

    % Extract center profile
    profile = results.frames(center, :, i);
    
    % Determine location of maximum
    [~, I] = max(profile);
    
    % Circshift to center maximum
    profile = circshift(profile, floor(length(profile) / 2) - I, 2);
    theta = circshift(results.itheta, floor(length(profile) / 2) - I, 2);
    frame = squeeze(circshift(results.frames(:, :, i), ...
        floor(length(profile) / 2) - I, 2));
    if exist('Event', 'file') == 2
        Event(sprintf('Frame %i circshifted to center on position %i', ...
            i, I));
    end
    
    % Re-determine location and value of maximum
    [C, I] = max(profile);
    if exist('Event', 'file') == 2
        Event(sprintf(['Radial entrance profile maximum identified as %g', ...
            ' at angle %0.3f deg'], C, theta(1, I)));
    end

    % Find highest lower index just below half maximum
    lI = find(profile(1:I) < C/2, 1, 'last');

    % Find lowest upper index just above half maximum
    uI = find(profile(I:end) < C/2, 1, 'first');
    
    % Verify edges were found
    if isempty(uI) || isempty(lI)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf('Field edges were not found for frame %i', i), ...
                'WARN');
        end
    
        % Continue to next profile
        continue;
        
    % Otherwise, verify edges are sufficiently far from array edges
    elseif lI-1 < 1 || lI+2 > length(profile) || I+uI-3 < 1 || ...
            I+uI > length(profile)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf(['Field edges are too close to profile', ...
                ' edge to compute for frame %i'], i), 'WARN');
        end
        
        % Continue to next profile
        continue;
    end
        
    % Interpolate to find lower half-maximum value
    l = interp1(profile(lI-1:lI+2), theta(1, lI-1:lI+2), C/2, 'spline');

    % Interpolate to find upper half-maximum value
    u = interp1(profile(I+uI-3:I+uI), theta(1, I+uI-3:I+uI), C/2, 'spline');

    % Compute angle as average of l and u thetas
    if (u < l); u = u + 360; end
    results.alpha(1,i) = mod((u + l) / 2 - rot, 360);
    if exist('Event', 'file') == 2
        Event(sprintf(['Frame %i entrance beam angle computed as ', ...
            '%0.3f deg'], i, results.alpha(1,i)));
    end

    % Interpolate longitudinal profile
    long = interp1(1:size(frame,2), frame(:,:)', (I + uI + lI) / 2);

    % Determine location and value of longitudinal maximum
    [C, I] = max(long);
    if exist('Event', 'file') == 2
        Event(sprintf(['Longitudinal entrance profile maximum ', ...
            'identified as %g at Y = %0.3f cm'], C, results.iY(I,1)));
    end

    % Find highest lower index just below half maximum
    lI = find(long(1:I) < C/2, 1, 'last');

    % Find lowest upper index just above half maximum
    uI = find(long(I:end) < C/2, 1, 'first');

    % Verify longitudinal edges were found
    if isempty(uI) || isempty(lI)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf('Field edges were not found for frame %i', i), ...
                'WARN');
        end

    % Otherwise, verify edges are sufficiently far from array edges
    elseif lI-1 < 1 || lI+2 > length(long) || I+uI-3 < 1 || ...
            I+uI > length(long)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf(['Field edges are too close to profile', ...
                ' edge to compute for frame %i'], i), 'WARN');
        end

    % Otherwise, continue longitudinal edge calculation
    else

        % Interpolate to find lower half-maximum value
        l = interp1(long(lI-1:lI+2), results.iY(lI-1:lI+2,1), C/2, ...
            'spline');

        % Interpolate to find upper half-maximum value
        u = interp1(long(I+uI-3:I+uI), results.iY(I+uI-3:I+uI,1), C/2, ...
            'spline');

        % Store field center
        results.beta(1,i) = (u+l)/2;
        if exist('Event', 'file') == 2
            Event(sprintf(['Entrance Y FWHM defined center ', ...
                'identified at %0.3f cm'], results.beta(1,i)));
        end
    end
  
    % Clear temporary variables
    clear C I l lI u uI long; 
    
    % Circshift to center the exit profile
    profile = circshift(profile, 180, 2);
    theta = circshift(theta, 180, 2);
    frame = squeeze(circshift(frame, 180, 2));
    if exist('Event', 'file') == 2
        Event(sprintf('Frame %i circshifted 180 deg to exit profile', i));
    end
        
    % Determine location and value of maximum
    [C, I] = max(profile(floor(size(profile, 2) * 1/4):...
        floor(size(profile, 2) * 3/4)));
    I = I + floor(size(profile, 2) / 4);
    if exist('Event', 'file') == 2
        Event(sprintf(['Radial exit profile maximum identified as %g at', ...
            ' angle %0.3f deg'], C, theta(1,I)));
    end
    
    % Find highest lower index just below half maximum
    lI = find(profile(1:I) < C/2, 1, 'last');

    % Find lowest upper index just above half maximum
    uI = find(profile(I:end) < C/2, 1, 'first');
    
    % Verify edges were found
    if isempty(uI) || isempty(lI)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf('Field edges were not found for frame %i', i), ...
                'WARN');
        end
    
        % Continue to next profile
        continue;
        
    % Otherwise, verify edges are sufficiently far from array edges
    elseif lI-1 < 1 || lI+2 > length(profile) || I+uI-3 < 1 || ...
            I+uI > length(profile)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf(['Field edges are too close to profile', ...
                ' edge to compute for frame %i'], i), 'WARN');
        end
        
        % Continue to next profile
        continue;
    end
        
    % Interpolate to find lower half-maximum value
    l = interp1(profile(lI-1:lI+2), theta(1,lI-1:lI+2), C/2, 'spline');

    % Interpolate to find upper half-maximum value
    u = interp1(profile(I+uI-3:I+uI), theta(1,I+uI-3:I+uI), C/2, 'spline');

    % Compute angle as average of l and u thetas
    if (u < l); u = u + 360; end
    results.alpha(2,i) = mod((u + l) / 2 - rot, 360);
    if exist('Event', 'file') == 2
        Event(sprintf(['Frame %i exit beam angle computed as ', ...
            '%0.3f deg'], i, results.alpha(2,i)));
    end

    % Interpolate longitudinal profile
    long = interp1(1:size(frame,2), frame(:,:)', (I + uI + lI) / 2);
    if exist('Event', 'file') == 2
        Event(['Longitudinal exit profile interpolated along ', ...
            'profile center']);
    end

    % Determine location and value of longitudinal maximum
    [C, I] = max(long);
    if exist('Event', 'file') == 2
        Event(sprintf(['Longitudinal exit profile maximum ', ...
            'identified as %g at Y = %0.3f cm'], C, results.iY(I,1)));
    end

    % Find highest lower index just below half maximum
    lI = find(long(1:I) < C/2, 1, 'last');

    % Find lowest upper index just above half maximum
    uI = find(long(I:end) < C/2, 1, 'first');

    % Verify longitudinal edges were found
    if isempty(uI) || isempty(lI)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf('Field edges were not found for frame %i', i), ...
                'WARN');
        end

    % Otherwise, verify edges are sufficiently far from array edges
    elseif lI-1 < 1 || lI+2 > length(long) || I+uI-3 < 1 || ...
            I+uI > length(long)

        % Log event
        if exist('Event', 'file') == 2
            Event(sprintf(['Field edges are too close to profile', ...
                ' edge to compute for frame %i'], i), 'WARN');
        end

    % Otherwise, continue longitudinal edge calculation
    else

        % Interpolate to find lower half-maximum value
        l = interp1(long(lI-1:lI+2), results.iY(lI-1:lI+2,1), C/2, ...
            'spline');

        % Interpolate to find upper half-maximum value
        u = interp1(long(I+uI-3:I+uI), results.iY(I+uI-3:I+uI,1), C/2, ...
            'spline');

        % Store field center
        results.beta(2,i) = (u+l)/2;
        if exist('Event', 'file') == 2
            Event(sprintf(['Exit Y FWHM defined center ', ...
                'identified at %0.3f cm'], results.beta(2,i)));
        end
    end
    
    % Clear temporary variables
    clear C I l lI u uI profile theta frame long;
end

% Log completion
if exist('Event', 'file') == 2
    
    % Log event 
    Event(sprintf(['SNC ArcCHECK analysis successfully completed in ', ...
      '%0.3f seconds'], toc(timer)));
  
    % Clear temporary variable
    clear timer;
end

% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

% Clear temporary variables
clear center c i progress rot;

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
