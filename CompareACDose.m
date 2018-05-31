function results = CompareACDose(varargin)
% AnalyzeACDose compares measured dose from ArcCHECK file exports (either
% .acm or .txt) to matching reference doses file and computes agreement 
% statistics, including gamma pass rate using the provided Gamma criteria. 
% The tool accepts two primary inputs; a folder containing the measurement
% files, and a second folder containing the reference DICOM RTPLAN and
% RTDOSE files. The tool will scan each folder, and attempt to match the
% DICOM RT patient name/ID and plan name to the measurement file name. If a
% match is found, the tool will run the comparison and store the results in
% the Plans field of the returned structure.
%
% As part of the analysis, this function accepts a Name/Value configuration
% option 'PlanClasses' with an n x 2 cell array of Plan class names and
% associated case insensitive regular expressions to match to the Plan 
% name. If provided, the function will attempt to categorize each 
% measurement into a Plan class, then compare each Plan type. Note, each 
% Plan is only matched to the first class that matches the regexp.
%
% The following required and optional inputs are available:
%   varargin{1}: file/folder/list of files containing measurement files
%       (see ParseD4tables for information on compatible formats)
%   varargin{2}: file/folder/list of files containing reference DICOM RT
%       Plan and RT Dose files for each measurement file in varargin{1}
%   varargin{3}: (optional) reserved for future use. For now pass anything.
%   varargin{4:nargin}: Name/Value pairs of analysis criteria, such as
%       reference dose mode 'RefDose' or Gamma settings ('GammaAbs', 
%       'GammaDTA','GammaRange'). See below for a full list of options. A
%       plan class cell array can also be provided using 'PlanClasses', 
%       as described above.
%
% Upon successful completion, a structure is returned with these fields:
%   GammaAbs: vector of absolute Gamma criteria evaluated, as a percent
%   GammaDTA: vector of absolute Gamma criteria evaluated (together, with
%       GammaAbs, this produces a 2D table of all permutations), in mm
%   GammaRange: 2 element vector of the low and high percentage of the
%       reference dose to be included in Gamma statistics
%   GammaLimit: the number of DTAs to search when computing Gamma (for
%       example, is GammaLimit is 2 and GammaDTA is 3, the code will search
%       in a 6 mm radius around each measurement point)
%   GammaRes: Gamma DTA search resolution (a value of 10 and GammaDTA of 3
%       mm means that the code will evaluate every 0.3 mm).
%   RefDose: a string indicating the type of reference dose. Can be
%       'measmax' (the reference dose will be the maximum measured dose) or
%       'planmax' (the reference dose will be the maximum reference dose)
%   AbsRange: 2 element vector of the low and high percentage of the
%       reference dose to be included in absolute dose difference
%   Plans: table with one row for each Plan matched to a dose volume
%   
% The following example illustrates how to use this function:
%
% % Define folder containing ArCHECK measurements
% meas = '~/Documents/QA Systems/ArcCHECK/Tomo Measurement Example/';
%
% % Define folder containing reference dose (can be the same folder)
% tps = '~/Documents/QA Systems/ArcCHECK/Tomo Measurement Example/';
%
% % Execute evaluation using default criteria
% results = CompareACDose(meas, tps)
%
% % Execute evaluation, using folders above and 2%/2mm criteria
% results = CompareACDose(meas, tps, [], 'GammaAbs', 2, 'GammaDTA', 2);
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2018 University of Wisconsin Board of Regents
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

% Set default options
results.Progress = true;
results.AbsRange = [50 500];
results.GammaAbs = 3;
results.GammaDTA = 3;
results.GammaRange = [20 500];
results.GammaLimit = 2;
results.GammaRes = 10;
results.RefDose = 'measmax';
results.PlanClasses = cell(0,2);

% Update options
for i = 4:2:nargin
    results.(varargin{i}) = varargin{i+1};
end

% Add dicom_tools submodule to search path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(path, 'dicom_tools'));

% Check if MATLAB can find ScanDICOMPath
if exist('ScanDICOMPath', 'file') ~= 2

    % If not, throw an error
    if exist('Event', 'file') == 2
        Event(['The dicom_tools submodule does not exist in the search path. ', ...
            'Use git clone --recursive or git submodule init followed by git ', ...
            'submodule update to fetch all submodules'], 'ERROR');
    else
        error(['The dicom_tools submodule does not exist in the search path. ', ...
            'Use git clone --recursive or git submodule init followed by git ', ...
            'submodule update to fetch all submodules']);
    end
end

%% Scan reference dose folders for DICOM files
% Display progress bar
if usejava('jvm') && feature('ShowFigureWindows') && results.Progress
    progress = waitbar(0, 'Scanning reference dataset for RT Plans');
end

% Scan first reference folder
if exist('Event', 'file') == 2
    Event('Scanning first reference dataset for RT Plans');
    t = tic;
end

% Store cell array of DICOM files as ref{1}
ref{1} = ScanDICOMPath(varargin{2}, 'Progress', false);

% Log start of measurement folders
if exist('Event', 'file') == 2
    Event('Scanning for measurement files');
end

% Scan measurements folder directory contents
if iscell(varargin{1})
    meas = varargin{1};
elseif isfolder(varargin{1})
    meas = dir(fullfile(varargin{1}, '**'));
else
    meas = dir(varargin{1});
end

% Generate IEC X/Y/Z coordinates
[z, x, y] = pol2cart(repmat(0:pi/33:(2*pi-pi/33), 1, 21) - pi, ...
    33/pi, -10.5:1/66:(10.5-1/66));

%% Initialize data structures
% Initialize Plan list, beam list, and gamma arrays return fields
fixedNames = {'MeasurementFile', 'Patient', 'ID', 'Plan', 'Machine', ...
    'Energy', 'PlanClass'};
fixedUnits = {'', '', '', '', '', 'MV', ''};
varNames = {'ReferenceFile', 'ReferenceDose', 'GammaPassRateGlobal', ...
    'GammaMeanGlobal', 'GammaMaxGlobal', 'GammaPassRateLocal', ...
    'GammaMeanLocal', 'GammaMaxLocal', 'MedianAbsDiff', 'AbsVolume', ...
    'GammaVolume'};
varUnits = {'', 'Gy', '%', '', '', '%', '', '', '%', 'cc', 'cc'};

% If only one reference exists
if length(ref) == 1
    v = varNames;
    u = varUnits;

% If multiple references exist
else
    
    % Initialize variable names
    v = cell(1,length(varNames)*length(ref));
    u = cell(1,length(varNames)*length(ref));
    for i = 1:length(ref)
        for j = 1:length(varNames)
            
            % Append a letter (A, B, C, ...) onto each variable name
            v{(i-1)*length(varNames)+j} = [varNames{j}, char(i+64)];
            u{(i-1)*length(varNames)+j} = varUnits{j};
        end
    end
end
results.Datasets = length(ref);

% Create tables
results.Plans = array2table(zeros(0,length(fixedNames)+length(v)), ...
    'VariableNames', horzcat(fixedNames, v));
results.Plans.Properties.VariableUnits = horzcat(fixedUnits, u);

% Clear and initialize GPU memory.  If CUDA is not enabled, or if the
% Parallel Computing Toolbox is not installed, this will error, and the
% function will automatically rever to CPU computation via the catch
% statement
try
    gpuDevice(1);
catch
    if exist('Event', 'file') == 2
        Event('GPU failed, will perform CPU interpolation', 'WARN');
    end
end

% Loop through each file
for i = 1:length(meas)
    
    % If the folder content is . or .. or a folder, skip to next file
    if isstruct(meas(i)) && (strcmp(meas(i).name, '.') || ...
            strcmp(meas(i).name, '..') || meas(i).isdir == 1)
        continue;
    end
    
    % Update waitbar
    if exist('progress', 'var') && ishandle(progress) && results.Progress
        waitbar(0.2 + 0.7*i/length(meas), progress, ...
            sprintf('Comparing ArcCheck measurement %i of %i', i, ...
            length(meas)));
    end

    % Parse file extension
    if isstruct(meas(i))
        [path, name, ext] = fileparts(fullfile(meas(i).folder, meas(i).name));
    elseif iscell(meas(i))
        [path, name, ext] = fileparts(meas{i});
    else
        [path, name, ext] = fileparts(meas(i));
    end
    
    %% Parse measurement
    try

        % Parse .acm file
        if strcmpi(ext, '.acm')
            m = ParseSNCacm(path, [name, ext], 'Progress', false);

            % Compute background and calibration-corrected dose
            data = table(m.x, m.y, m.z, m.theta, ...
                (m.data(end, 12:1397) - m.data(end,3) .* ...
                m.background(2:1387)) .* m.calibration(2:1387) * ...
                m.dosecal/100, 'VariableNames', {'x', 'y', 'z', 'theta', ...
                'measured'});

        % Parse .txt file
        elseif strcmpi(ext, '.txt')
            m = ParseSNCtxt(path, [name, ext], 'Type', 'ARCCHECK', ...
                'Progress', false);

            % Use dose counts for measured dose
            dc = flip(m.dosecounts',2);
            data = table(x', y', z', atan2d(z, x)' + 90, ...
                dc(dc > 0)/100, 'VariableNames', {'x', 'y', 'z', 'theta', ...
                'measured'});
            clear dc;
            
        else
            continue;
        end
        
    % If parser fails, catch gracefully
    catch
        if exist('Event', 'file') == 2
            Event('Parsing failed, skipping measurement', 'WARN');
        else
            warning(['Parsing failed, skipping ', name, ext]);
        end
        continue;
    end
    
    % Initialize plan IDs
    matched = zeros(1, length(ref));
    
    %% Match plan to reference dose files
    % Loop through each reference cell array
    for r = 1:length(ref)
        for j = 1:length(ref{r})
        
            % Split name
            nr = strsplit(ref{r}{j,5}, ',');
            
            % If reference is a plan, and both the plan name and either the
            % patient name or ID matches the file name (plan information is
            % not stored in the SNC exported files)
            if strcmp(ref{r}{j,10}, 'PLAN') && ...
                    contains(name, ref{r}{j,9}, 'IgnoreCase', true) && ...
                    (contains(name, ref{r}{j,6}, 'IgnoreCase', true) || ...
                    (contains(name, nr{1}, 'IgnoreCase', true) && ...
                    (length(nr) < 2 || contains(name, nr{2}, ...
                    'IgnoreCase', true))))
                matched(r) = j;
                break;
            end
        end
    end
    
    % If all references were matched
    if all(matched)
        
        % If all beams have the same energy, use that, otherwise report Mixed
        if length(unique(ref{1}{matched(1),13})) == 1
            energy = ref{1}{matched(1),13}{1};
        else
            energy = 'Mixed';
        end
        
        % Match plan class (will stop after the first match is found)
        class = 'Undefined';
        for j = 1:size(results.PlanClasses, 1)
            if ~isempty(regexpi(ref{1}{matched(1),9}, ...
                    results.PlanClasses{j,2}))
                class = results.PlanClasses{j,1};
                break;
            end
        end
        
        % Initialize results
        result = {[name, ext], ref{1}{matched(1),5}, ref{1}{matched(1),6}, ...
            ref{1}{matched(1),9}, ref{r}{matched(1),12}{1}, energy, class};
        
        % Loop through each reference
        for j = 1:length(matched)
            
            % Load dose
            dose = LoadDICOMDose(ref{j}{matched(j),2}, ref{j}{matched(j),1});
            
            % Set reference value
            switch results.RefDose
                case 'measmax'
                    refval = max(data.measured);
                case 'planmax'
                    refval = max(max(max(dose.data)));
                otherwise
                    refval = 0;
            end
            
            % Compute dose volumes
            absvol = (sum(sum(sum(dose.data >= refval * ...
                results.AbsRange(1)/100))) - sum(sum(sum(dose.data > ...
                refval * results.AbsRange(2)/100)))) * prod(dose.width);
            gamvol = (sum(sum(sum(dose.data >= refval * ...
                results.GammaRange(1)/100))) - sum(sum(sum(dose.data > ...
                refval * results.GammaRange(2)/100)))) * prod(dose.width);
            
            % Log calculation
            if exist('Event', 'file') == 2
                Event(['Computing plan Gamma table between ', name, ext, ...
                    'and ', ref{j}{matched(j),1}]);
            end
            
            % Apply red laser shift to dose coordinates
            if ~isempty(ref{j}{matched(j),14}) && ...
                    length(ref{j}{matched(j),14}) == 3
                dose.start = dose.start - [ref{j}{matched(j),14}(1), ...
                    ref{j}{matched(j),14}(3), -ref{j}{matched(j),14}(2)]/10;
            end
            
            % Append reference file and gamma table
            result = [result, ref{j}{matched(j),1}, refval, ...
                gammaTable(data, dose, refval, results), absvol, ...
                gamvol]; %#ok<*AGROW>
        end
        
        % Append result onto results
        results.Plans = [results.Plans; result];
    
    % Otherwise, warn the user
    else
        if exist('Event', 'file') == 2
            Event(['No matching reference files were found for ', ...
                name, ext], 'WARN');
        else
            warning(['No matching reference files were found for ', ...
                name, ext])
        end
    end
end
    
%% Finish up
% Close waitbar
if exist('progress', 'var') == 1 && ishandle(progress)
    waitbar(1.0, progress);
    close(progress);
end

% Remove execution flags from results
results = rmfield(results, 'Progress');

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Comparison successfully completed in %0.3f seconds, ', ...
        'scanning %i Plans'], toc(t), size(results.Plans,1)));
end

% Clear temporary variables
clear t;

end

%% Compute Gamma Table Subfunction
function stats = gammaTable(meas, dose, refval, criteria)
% gammaTable is called by CompareD4Dose and computes a 2D gamma table based
% on a provided measurement array and reference dose volume, using the
% results GammaAbs, GammaDTA, GammaRange, GammaRes, and GammaLimit fields.

% Compute mesh grid for reference dose
[meshA, meshB, meshC] = meshgrid(single(dose.start(2) + ...
    dose.width(2) * (size(dose.data,2) - 1):-dose.width(2):dose.start(2)), ...
    single(dose.start(1):dose.width(1):dose.start(1) + dose.width(1)...
    * (size(dose.data,1) - 1)), single(dose.start(3):dose.width(3):...
    dose.start(3) + dose.width(3) * (size(dose.data,3) - 1)));

% Interpolate dose at original position
try
    % Run GPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = gather(interp3(gpuArray(meshA), gpuArray(meshB), ...
        gpuArray(meshC), gpuArray(single(dose.data)), ...
        gpuArray(meas.z), gpuArray(meas.x), gpuArray(-meas.y), ...
        'linear', 0));
catch

    % Run CPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = interp3(meshA, meshB, meshC, single(dose.data), meas.z, meas.x, ...
        -meas.y, '*linear', 0);
end

% Compute the median absolute dose difference
absDiff = (meas.measured - i)./refval * 100;
stats = cell(1,7);
stats{7} = median(absDiff(meas.measured > refval * criteria.AbsRange(1)/100 & ...
    meas.measured < refval * criteria.AbsRange(2)/100));

% Apply gamma range to measured values
meas = meas(meas.measured > refval * criteria.GammaRange(1)/100 & ...
    meas.measured < refval * criteria.GammaRange(2)/100,:);

% Compute shifts (in cm)
shifts = repmat(-criteria.GammaRes * criteria.GammaLimit:criteria.GammaRes * ...
    criteria.GammaLimit, size(meas,1), 1) / ...
    criteria.GammaRes * max(criteria.GammaDTA)/10;

% Compute shifted 2D arrays for each position
measA = single(horzcat(repmat(meas.z, 1, size(shifts,2)) + shifts, ...
    repmat(meas.z, 1, 2*size(shifts,2))));
measB = single(horzcat(repmat(meas.x, 1, size(shifts,2)), ...
    repmat(meas.x,1, size(shifts,2)) + shifts, ...
    repmat(meas.x, 1, size(shifts,2))));
measC = single(horzcat(repmat(-meas.y, 1, 2*size(shifts,2)), ...
    repmat(-meas.y, 1, size(shifts,2)) + shifts));

% Interpolate dose at shifted position
try
    % Run GPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = gather(interp3(gpuArray(meshA), gpuArray(meshB), ...
        gpuArray(meshC), gpuArray(single(dose.data)), ...
        gpuArray(measA), gpuArray(measB), gpuArray(measC), 'linear', 0));
catch

    % Run CPU interp3 function to compute the dose values at the shifted 
    % measured coordinate points
    i = interp3(meshA, meshB, meshC, single(dose.data), measA, measB, ...
        measC, '*linear', 0);
end

% Compute local gamma as min of each shifted position
loc = zeros(size(i,1), length(criteria.GammaAbs), length(criteria.GammaDTA));
glob = zeros(size(i,1), length(criteria.GammaAbs), length(criteria.GammaDTA));

% Loop through Abs, DTA criteria
for a = 1:length(criteria.GammaAbs)
    for d = 1:length(criteria.GammaDTA)
        
        % Compute global gamma based on refval as min of each shifted position
        glob(:,a,d) = sqrt(min(((repmat(meas.measured,1,size(i,2)) - i) ./ ...
            (refval * criteria.GammaAbs(a)/100)).^2 + ...
            (repmat(shifts,1,3)/(criteria.GammaDTA(d)/10)).^2,[],2));
        
        % Compute local gamma as min of each shifted position
        loc(:,a,d) = sqrt(min(((repmat(meas.measured,1,size(i,2)) - i) ./ ...
            (i * criteria.GammaAbs(a)/100)).^2 + ...
            (repmat(shifts,1,3)/(criteria.GammaDTA(d)/10)).^2,[],2));
    end
end

% Compute pass rate, mean, and maximum Gamma statistics for global gamma
stats{1} = zeros(length(criteria.GammaAbs), length(criteria.GammaDTA));
for a = 1:length(criteria.GammaAbs)
    for d = 1:length(criteria.GammaDTA)
        stats{1}(a,d) = sum(glob(:,a,d) <= 1) / size(glob,1) * 100;
    end
end
stats{2} = squeeze(mean(glob,1));
stats{3} = squeeze(max(glob,[], 1));

% Compute pass rate, mean, and maximum Gamma statistics for local gamma
stats{4} = zeros(length(criteria.GammaAbs), length(criteria.GammaDTA));
for a = 1:length(criteria.GammaAbs)
    for d = 1:length(criteria.GammaDTA)
        stats{4}(a,d) = sum(loc(:,a,d) <= 1) / size(loc,1) * 100;
    end
end
stats{5} = squeeze(mean(loc,1));
stats{6} = squeeze(max(loc,[], 1));

% Clear temporary variables
clear measA measB measC a d i shifts glob loc meshA meshB meshC absDiff;

end