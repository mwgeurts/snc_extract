function data = ParseSNCxml(path, names)
% ParseSNCxml reads in SNC 3D SCANNER Water Tank XML exported data files
% and creates a MATLAB structure of the results. The scan data is stored
% under the structure field data.Scans{i}.Layers{j}.Readings, with columns
% for each numerical value of the <ScanReading> tag (Reading ID, Sequence,
% X, Y, Z, Depth, Seconds, Charge, Current, and RelativeDose). If multiple
% files are provided, the Scans array is concatenated sequentially.
%
% Generally speaking, each MATLAB structure field follows the name/value of
% each <p> tag in the file. In order to improve efficiency, I have removed
% some of the nested containers by storing their values in the higher level
% structure (this reduces the number of levels in the structure and makes
% it overall easier to navigate and find things). I've also made certain
% assumptions about the order of certain tag groups (X/Y/Z/etc.) in order
% to improve the efficiency of this function.
%
% The input variable path should be a string containing the relative or
% absolute path to each SNC file, while names can either be a string
% containing the file name or a cell array of strings for each file name.
%
% Upon successful completion, below is a high level summary of the returned
% MATLAB structure (note not all fields are shown). The function attempts
% to identify the format of each returned field as either numerical, text,
% logical, or date/time (stored as a MATLAB datenum).
%
% data.RadiationDevices (cell array)
%   data.RadiationDevices.Energies (cell array)
%   data.RadiationDevices.Fields (cell array)
% data.Scans (cell array)
%   data.Scans.Layers (cell array)
%       data.Scans.Layers.Readings (table)
%   data.Scans.FieldSize (structure)
%   data.Scans.Electrometer (structure)
%   data.Scans.Detectors (structure)
%       data.Scans.Detectors.Field (structure)
%       data.Scans.Detectors.Reference (structure)
%   data.Scans.Background (structure)
%       data.Scans.Detectors.Field (structure)
%       data.Scans.Detectors.Reference (structure)
%   data.Scans.Collimation (structure)
%   data.Scans.Points (structure)
%   data.Scans.MeasurementMode (structure)
% data.version (double)
% data.exported (MATLAB datenum)
% data.nodes (cell array of strings)
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2017 University of Wisconsin Board of Regents
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
    Event(['Loading SNC water tank file ', strjoin(names, ...
        '\nLoading SNC water tank file ')]);
    tic;
end

% If a valid screen size is returned (MATLAB was run without -nodisplay)
if usejava('jvm') && feature('ShowFigureWindows')
    
    % Start waitbar
    progress = waitbar(0, 'Parsing SNC Water Tank file(s)');
end

% Initialize return structure
data.RadiationDevices = cell(0);
data.Scans = cell(0);

% The patient XML is parsed using xpath class
import javax.xml.xpath.*

% Loop through each file
for f = 1:length(names)
    
    % Read in the patient XML and store the Document Object Model node
    if exist('Event', 'file') == 2
        Event('Loading file contents data using xmlread');
    end
    doc = xmlread(fullfile(path, names{f}));

    % Initialize a new xpath instance to the variable factory
    factory = XPathFactory.newInstance;

    % Initialize a new xpath to the variable xpath
    xpath = factory.newXPath;

    % Search for file export version and date
    expression = xpath.compile('snc_fileexport');

    % Evaluate xpath expression and retrieve the results
    nodeList = expression.evaluate(doc, XPathConstants.NODESET);

    % If snc_fileexport tag was found and contains a version attribute
    if nodeList.getLength > 0 && nodeList.item(0).hasAttribute('version')

        % Store version
        data.version(f) = ...
            str2double(nodeList.item(0).getAttribute('version'));
    else

        % Otherwise, warn the user that patient info wasn't found
        if exist('Event', 'file') == 2
            Event(['File export version could not be found. Verify that ', ...
                'this file is an SNC water tank XML file.'], 'ERROR');
        else
            error(['File export version could not be found. Verify that ', ...
                'this file is an SNC water tank XML file.']);
        end
    end

    % If an export date attribute exists
    if nodeList.item(0).hasAttribute('date')
        try
            data.exported(f) = ...
                datenum(char(nodeList.item(0).getAttribute('date')), ...
                'mm/dd/yyyy HH:MM:SS');
        catch
            data.exported(f) = [];
        end
    end
    
    % Search for notes
    expression = xpath.compile('//snc_fileexport/notes');

    % Evaluate xpath expression and retrieve the results
    nodeList = expression.evaluate(doc, XPathConstants.NODESET);
    
    % If notes exist
    if nodeList.getLength > 0
        
        % Store notes
        data.notes{f} = char(nodeList.item(0).getFirstChild.getNodeValue);
    end
    
    % Search for radiation device information
    expression = xpath.compile(['//snc_fileexport/entityinfo/device/', ...
         'RadiationDevice']);

    % Evaluate xpath expression and retrieve the results
    nodeList = expression.evaluate(doc, XPathConstants.NODESET);
    
    % If device information exists
    for i = 1:nodeList.getLength

        % Retrieve nodes
        nodes = nodeList.item(i-1).getChildNodes;

        % Initialize temporary structure for RadiationDevice
        r = struct;
        
        % Loop through nodes
        for j = 0:nodes.getLength-1

            % Parse name and value
            if nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    nodes.item(j).hasAttribute('value')

                % Store parsed value
                r.(char(nodes.item(j).getAttribute('name'))) = ...
                    parse(char(nodes.item(j).getAttribute('value')));
                

            % Otherwise, if this is a collection of energies
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && strcmp('Energies', ...
                    char(nodes.item(j).getAttribute('name')))

                % Search for energies
                subexpression = xpath.compile(...
                    'RadiationDeviceEnergyCollection/RadiationDeviceEnergy');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Initialize energy cell array
                r.Energies = cell(subnodeList.getLength, 1);

                % Loop through energies
                for k = 1:subnodeList.getLength

                    % Initialize struct
                    r.Energies{k} = struct;

                    % Retrieve subnodes
                    subnodes = subnodeList.item(k-1).getChildNodes;

                    % Loop through subnodes
                    for l = 0:subnodes.getLength-1

                        % Parse name and value
                        if subnodes.item(l).hasAttributes && ...
                                subnodes.item(l).hasAttribute('name') && ...
                                subnodes.item(l).hasAttribute('value')
                           
                            % Store parsed value
                            r.Energies{k}.(char(subnodes.item(l)...
                                .getAttribute('name'))) = parse(char(...
                                subnodes.item(l).getAttribute('value')));
                        end
                    end
                end

            % Otherwise, if this is a collection of fields
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && strcmp('Fields', ...
                    char(nodes.item(j).getAttribute('name')))

                % Search for fields
                subexpression = xpath.compile('FieldCollection/Field');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Initialize energy cell array
                r.Fields = cell(subnodeList.getLength, 1);

                % Loop through energies
                for k = 1:subnodeList.getLength

                    % Initialize struct
                    r.Fields{k} = struct;

                    % Retrieve subnodes
                    subnodes = subnodeList.item(k-1).getChildNodes;

                    % Loop through subnodes
                    for l = 0:subnodes.getLength-1

                        % Parse name and value
                        if subnodes.item(l).hasAttributes && ...
                                subnodes.item(l).hasAttribute('name') && ...
                                subnodes.item(l).hasAttribute('value')
                            
                            % Store parsed value
                            r.Fields{k}.(char(subnodes.item(l)...
                                .getAttribute('name'))) = parse(char(...
                                subnodes.item(l).getAttribute('value')));
                        end
                    end
                end
                
            % Otherwise, if this is the institution
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp('Institution', ...
                    char(nodes.item(j).getAttribute('name')))
                
                % Search for institution information
                subexpression = xpath.compile('Institution/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);
                
                % Initialize institution
                r.Institution = struct;
                
                % Loop through fields
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        r.Institution.(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
                
                % Search for address information
                subexpression = xpath.compile('Institution/p/Location/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through fields
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        r.Institution.(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
                
                % Search for country information
                subexpression = ...
                    xpath.compile('Institution/p/Location/p/Country/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);
                
                % Loop through fields
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            strcmp(char(subnodeList.item(k)...
                            .getAttribute('name')), 'Name')
                       r.Institution.Country = ...
                            char(subnodeList.item(k).getAttribute('value'));
                    end
                end
            end
        end
        
        % Look for existing RadiationDevices with same UniqueID
        e = false;
        for j = 1:length(data.RadiationDevices)
            
            % If the UniqueID matches an existing device, flag it
            if isfield(data.RadiationDevices{j}, 'UniqueId') && ...
                    isfield(r, 'UniqueId') && ...
                    strcmp(data.RadiationDevices{j}.UniqueId, r.UniqueId)
                e = true;
                break;
            end
        end
        
        % If no matching UniqueIDs were gound
        if ~e
            data.RadiationDevices{length(data.RadiationDevices)+1} = r;
        end
        
        % Clear temporary variables
        clear e j k l r nodes subexpression subnodeList subnodes;
    end
    
    % Store index of existing scans
    num = length(data.Scans);
    
    % Search for scan entity information
    expression = xpath.compile('snc_fileexport/scans/scan/entityinfo');

    % Evaluate xpath expression and retrieve the results
    nodeList = expression.evaluate(doc, XPathConstants.NODESET);
    
    % Loop through scan entities
    for i = 1:nodeList.getLength
        
        % Store the institution
        if nodeList.item(i-1).hasAttribute('Institution')
            data.Scans{num+i}.Institution = ...
                char(nodeList.item(i-1).getAttribute('Institution'));
        end
        
        % Store the radiation device
        if nodeList.item(i-1).hasAttribute('RadiationDevice')
            data.Scans{num+i}.RadiationDevice = ...
                char(nodeList.item(i-1).getAttribute('RadiationDevice'));
        end
    end
    
    % Search for scan data
    expression = xpath.compile('snc_fileexport/scans/scan/Scan');

    % Evaluate xpath expression and retrieve the results
    nodeList = expression.evaluate(doc, XPathConstants.NODESET);
    
    % Loop through scan data
    for i = 1:nodeList.getLength
        
        % Update waitbar
        if exist('progress', 'var') && ishandle(progress)
            waitbar((f-1+(i-1)/nodeList.getLength)/length(names), progress);
        end
        
        % Retrieve nodes
        nodes = nodeList.item(i-1).getChildNodes;

        % Loop through nodes
        for j = 0:nodes.getLength-1
            
            % If node contains name/value information
            if nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    nodes.item(j).hasAttribute('value')
                
                % Store parsed value
                data.Scans{num+i}.(char(nodes.item(j)...
                    .getAttribute('name'))) = parse(char(nodes.item(j)...
                    .getAttribute('value')));
                
            % Otherwise, if this contains Layer data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Layers')
                
                % Search for layers
                subexpression = ...
                    xpath.compile('ScanLayerCollection/ScanLayer');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);
                
                % Initialize layers with result
                data.Scans{num+i}.Layers = cell(subnodeList.getLength, 1);
                
                % Loop through layers
                for k = 1:subnodeList.getLength
                    
                    % Get layer subnodes
                    subnodes = subnodeList.item(k-1).getChildNodes;
                    
                    % Loop through subnodes
                    for l = 0:subnodes.getLength-1
                        
                        % Parse name and value
                        if subnodes.item(l).hasAttributes && ...
                                subnodes.item(l).hasAttribute('name') && ...
                                subnodes.item(l).hasAttribute('value')
                            
                            % Store parsed value
                            data.Scans{num+i}.Layers{k}.(char(...
                                subnodes.item(l).getAttribute('name'))) = ...
                                parse(char(subnodes.item(l)...
                                .getAttribute('value')));
                        end
                    end
                    
                    % Search for Readings
                    subexpression = ...
                        xpath.compile('p/ScanReadingCollection/ScanReading');

                    % Evaluate xpath expression and retrieve the results
                    readings = subexpression.evaluate(subnodeList.item(k-1), ...
                        XPathConstants.NODESET);
                
                    % Initialize temporary readings array
                    arr = zeros(readings.getLength, 13);
                    
                    % Loop through subnodes
                    for l = 1:readings.getLength
                        
                        % Get reading subnodes
                        subnodes = readings.item(l-1).getChildNodes;
                        
                        % Initialize counter
                        c = 1;
                        
                        % Loop through subnodes
                        for m = 0:subnodes.getLength-1

                            % Store subnode data, assuming the order. Note,
                            % this is a little risky but significantly
                            % reduces the load time
                            if subnodes.item(m).hasAttributes
                                
                                % If this is the third child, these
                                % are textual data
                                if c == 3
                                    
                                    % Increment counter
                                    c = c + 1;
                                    
                                    % Skip these values
                                    continue;
                                
                                % If this is the fourth child, assume it is
                                % point data
                                elseif c == 4
                                    
                                    % Store DataPoint children
                                    p = subnodes.item(m).getChildNodes;

                                    % Loop through child nodes
                                    for n = 0:p.getLength-1

                                        % If this is the DataPoint child
                                        if p.item(n).hasAttributes
                                            p = p.item(n).getChildNodes;
                                            break;
                                        end
                                    end

                                    % Loop through DataPoint's children
                                    for n = 0:p.getLength-1
                                        
                                        % If this child has attributes
                                        if p.item(n).hasAttributes
                                           
                                            % Store value
                                            arr(l,c) = str2double(p.item(n)...
                                                .getAttribute('value'));
                                            
                                            % Increment counter
                                            c = c + 1;
                                        end
                                    end
                                
                                % Otherwise, if the 11th, stop searching
                                elseif c == 12
                                    break;
                                    
                                % Otherwise, store its numerical value
                                else
                                    
                                    % Store value
                                    arr(l,c) = str2double(subnodes.item(m)...
                                        .getAttribute('value'));
                                
                                    % Increment counter
                                    c = c + 1;
                                end
                            end
                        end
                    end
                    
                    % Convert array to table
                    data.Scans{num+i}.Layers{k}.Readings = array2table(arr, ...
                        'VariableNames', {'ReadingId', 'Sequence', 'Layer', ...
                        'X', 'Y', 'Z', 'Depth', 'Seconds', 'Charge', 'Current', ...
                        'RelativeDose', 'Name', 'Modified'});
                    
                    % Clear temporary variables
                    clear l c n m p arr subexpression readings subnodes;
                end
                
                % Clear temporary variables
                clear k subexpression subnodeList;
                
            % Otherwise, if this contains field size data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'FieldSize')
                
                % Search for field data
                subexpression = xpath.compile('Field/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through field data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.FieldSize.(char(subnodeList...
                            .item(k).getAttribute('name'))) = ...
                            parse(char(subnodeList.item(k)...
                            .getAttribute('value')));
                    end
                end
                
                % Clear temporary variables
                clear k subexpression subnodeList;
                
            % Otherwise, if this contains electrometer data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Electrometer')
                
                % Search for electrometer data
                subexpression = xpath.compile('ScanElectrometer/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through field data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.Electrometer.(char(subnodeList...
                            .item(k).getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                        
                    % Parse container
                    elseif subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name')
                        name = ...
                            char(subnodeList.item(k).getAttribute('name'));
                        
                        % Store temperatures
                        if contains(name, 'Temperature')
                        
                            % Search for electrometer data
                            subexpression = xpath.compile('Temperature/p');

                            % Evaluate xpath expression and retrieve the results
                            t = subexpression.evaluate(subnodeList.item(k), ...
                                XPathConstants.NODESET);
                            
                            % Store Amount
                            for l = 0:t.getLength-1
                                if strcmp(char(t.item(l)...
                                        .getAttribute('name')), 'Amount')
                                    data.Scans{num+i}.Electrometer.(name) = ...
                                        str2double(t.item(l)...
                                        .getAttribute('value'));
                                end
                            end
                            
                        % Store pressures
                        elseif contains(name, 'Pressure')
                            
                            % Search for electrometer data
                            subexpression = xpath.compile('Pressure/p');

                            % Evaluate xpath expression and retrieve the results
                            t = subexpression.evaluate(subnodeList.item(k), ...
                                XPathConstants.NODESET);
                            
                            % Store Amount
                            for l = 0:t.getLength-1
                                if strcmp(char(t.item(l)...
                                        .getAttribute('name')), 'Amount')
                                    data.Scans{num+i}.Electrometer.(name) = ...
                                        str2double(char(t.item(l)...
                                        .getAttribute('value')));
                                end
                            end
                        end
                    end
                end
                
                % Clear temporary variables
                clear k t l name subexpression subnodeList;
                
            % Otherwise, if this contains detector data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Detectors')
                
                % Search for detectors
                subexpression = xpath.compile(...
                    'ScanDetectorCollection/ScanDetector');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through detector data
                for k = 0:subnodeList.getLength-1
                    
                    % Search for type data
                    subexpression = xpath.compile('p');

                    % Evaluate xpath expression and retrieve the results
                    detector = subexpression.evaluate(...
                        subnodeList.item(k), ...
                        XPathConstants.NODESET);
                    
                    % Store detector type
                    for l = 0:detector.getLength-1
                        if detector.item(l).hasAttribute('name') && ...
                                strcmp(detector.item(l).getAttribute('name'), ...
                                'Type')
                            type = char(detector.item(l)...
                                .getAttribute('value'));
                        end
                    end
                    
                    % Store detector scan parameters
                    for l = 0:detector.getLength-1
                        if detector.item(l).hasAttribute('name') && ...
                                detector.item(l).hasAttribute('value')
                            
                            % Store parsed value
                            data.Scans{num+i}.Detectors.(type).(char(...
                                detector.item(l).getAttribute('name'))) = ...
                                parse(char(detector.item(l)...
                                .getAttribute('value')));
                        end
                    end
                    
                    % Search for Detector data
                    subexpression = xpath.compile('p/Detector/p');

                    % Evaluate xpath expression and retrieve the results
                    detector = subexpression.evaluate(...
                        subnodeList.item(k), ...
                        XPathConstants.NODESET);
                    
                    % Store detector data
                    for l = 0:detector.getLength-1
                        if detector.item(l).hasAttribute('name') && ...
                                detector.item(l).hasAttribute('value')
                            
                            % Store parsed value
                            data.Scans{num+i}.Detectors.(type).(char(...
                                detector.item(l).getAttribute('name'))) = ...
                                parse(char(detector.item(l)...
                                .getAttribute('value')));
                        end
                    end
                    
                    % Search for offset data
                    subexpression = xpath.compile('p/DataPointOffset/p');

                    % Evaluate xpath expression and retrieve the results
                    detector = subexpression.evaluate(...
                        subnodeList.item(k), ...
                        XPathConstants.NODESET);
                    
                    % Store detector data
                    arr = zeros(detector.getLength,1);
                    for l = 1:detector.getLength
                        arr(l) = str2double(detector.item(l-1)...
                            .getAttribute('value'));
                    end
                    data.Scans{num+i}.Detectors.(type).PointOffset = ...
                        array2table(arr', 'VariableNames', {'X', 'Y', 'Z'});
                end
                
                % Clear temporary variables
                clear k l detector type subexpression subnodeList arr;
                
            % Otherwise, if this contains background data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Background')
                
                % Search for detectors
                subexpression = xpath.compile(...
                    'ScanBackgroundCollection/ScanBackground');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through detector data
                for k = 0:subnodeList.getLength-1
                    
                    % Search for type data
                    subexpression = xpath.compile('p');

                    % Evaluate xpath expression and retrieve the results
                    detector = subexpression.evaluate(...
                        subnodeList.item(k), ...
                        XPathConstants.NODESET);
                    
                    % Store detector type
                    for l = 0:detector.getLength-1
                        if detector.item(l).hasAttribute('name') && ...
                                strcmp(detector.item(l).getAttribute('name'), ...
                                'Type')
                            type = char(detector.item(l)...
                                .getAttribute('value'));
                        end
                    end
                    
                    % Store detector background parameters
                    for l = 0:detector.getLength-1
                        if detector.item(l).hasAttribute('name') && ...
                                detector.item(l).hasAttribute('value')
                            
                            % Store parsed value
                            data.Scans{num+i}.Background.(type).(char(...
                                detector.item(l).getAttribute('name'))) = ...
                                parse(char(detector.item(l)...
                                .getAttribute('value')));
                        end
                    end
                    
                    % Search for Detector data
                    subexpression = xpath.compile('p/Detector/p');

                    % Evaluate xpath expression and retrieve the results
                    detector = subexpression.evaluate(...
                        subnodeList.item(k), ...
                        XPathConstants.NODESET);
                    
                    % Store detector background data
                    for l = 0:detector.getLength-1
                        if detector.item(l).hasAttribute('name') && ...
                                detector.item(l).hasAttribute('value')
                            
                            % Store parsed value
                            data.Scans{num+i}.Background.(type).(char(...
                                detector.item(l).getAttribute('name'))) = ...
                                parse(char(detector.item(l)...
                                .getAttribute('value')));
                        end
                    end
                end
                
                % Clear temporary variables
                clear k l detector type subexpression subnodeList;
                
            % Otherwise, if this contains wppb data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'WaterProofProfilerBackground')
                
                % Search for WPPB data
                subexpression = xpath.compile('WppBackground/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through WPPB data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.WaterProofProfilerBackground...
                            .(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
                
            % Otherwise, if this contains collimation data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Collimation')
                
                % Search for Collimation data
                subexpression = xpath.compile('ScanCollimation/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through Collimation data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.Collimation...
                            .(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
                
            % Otherwise, if this contains points data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Points')
                
                % Search for ScanDataPoint sets
                subexpression = xpath.compile('ScanDataPoint/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through point sets
                for k = 0:subnodeList.getLength-1

                    % If this is a name/value pair
                    if subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.Points.(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                        
                    % Otherwise, if a container
                    else
                        
                        % Search for DataPoint data
                        subexpression = xpath.compile('DataPoint/p');

                        % Evaluate xpath expression and retrieve the results
                        point = subexpression.evaluate(...
                            subnodeList.item(k), XPathConstants.NODESET);

                        % Store point data
                        arr = zeros(point.getLength, 1);
                        for l = 0:point.getLength-1
                            if point.item(l).hasAttribute('name') && ...
                                    point.item(l).hasAttribute('value')

                                % Store parsed value
                                arr(l+1) = str2double(char(point.item(l)...
                                    .getAttribute('value')));
                            end
                        end
                    
                        % Convert to table
                        data.Scans{num+i}.Points.(char(subnodeList.item(k)...
                            .getAttribute('name'))) = array2table(arr', ...
                            'VariableNames', {'X', 'Y', 'Z', 'Depth'});
                    end
                end
                
                % Clear temporary variables
                clear k l point arr subexpression subnodeList;
                
            % Otherwise, if this contains tpr data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'TissuePhantomRatio')
                
                % Search for TPR data
                subexpression = xpath.compile('ScanTissuePhantomRatio/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through TPR data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.TissuePhantomRatio...
                            .(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
                
            % Otherwise, if this contains wedge data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'Wedge')
                
                % Search for Wedge data
                subexpression = xpath.compile('ScanWedge/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through Wedge data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.Wedge...
                            .(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
                
            % Otherwise, if this contains mode data
            elseif nodes.item(j).hasAttributes && ...
                    nodes.item(j).hasAttribute('name') && ...
                    strcmp(char(nodes.item(j).getAttribute('name')), ...
                    'MeasurementMode')
                
                % Search for MeasurementMode data
                subexpression = xpath.compile('ScanMeasurementMode/p');

                % Evaluate xpath expression and retrieve the results
                subnodeList = subexpression.evaluate(nodes.item(j), ...
                    XPathConstants.NODESET);

                % Loop through MeasurementMode data
                for k = 0:subnodeList.getLength-1
                    
                    % Parse name and value
                    if subnodeList.item(k).hasAttributes && ...
                            subnodeList.item(k).hasAttribute('name') && ...
                            subnodeList.item(k).hasAttribute('value')
                        
                        % Store parsed value
                        data.Scans{num+i}.MeasurementMode...
                            .(char(subnodeList.item(k)...
                            .getAttribute('name'))) = parse(char(...
                            subnodeList.item(k).getAttribute('value')));
                    end
                end
            end
        end
        
        % Clear temporary variables
        clear j nodes
    end
    
    % Clear temporary variables
    clear i num expression nodeList doc factory xpath;
end

% Clear temporary variables
clear f;

% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
    clear progress;
end

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['%i radiation devices and %i scans successfully parsed ', ...
        'in %0.3f seconds'], length(data.RadiationDevices), ...
        length(data.Scans), toc));
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parsed = parse(value)
% parse is a subfunction of ParseSNCxml and attempts to parse a char array
% based on its data type. This function can find numerical values,
% logicals, and UTC date/time (ISO 8601)

% Parse empty/text values
if isempty(value)
    parsed = value;
    
% Parse numerical values
elseif regexp(value, '^-?\d*\.?\d*$', 'once')
    parsed = str2double(value);

% Parse logical values
elseif contains(value, 'true', 'IgnoreCase', true)
    parsed = true;

elseif contains(value, 'false', 'IgnoreCase', true)
    parsed = false;
    
% Parse UTC date/time values with microseconds
elseif regexp(value, '^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z$')
    try
        parsed = datenum(value, 'yyyy-mm-ddTHH:MM:SS.FFFZ');
    catch
        parsed = [];
    end

% Parse UTC date/time values without microseconds
elseif regexp(value, '^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\Z$')
    try
        parsed = datenum(value, 'yyyy-mm-ddTHH:MM:SSZ');
    catch
        parsed = [];
    end
    
% Parse textual values    
else
    parsed = value;
end

end