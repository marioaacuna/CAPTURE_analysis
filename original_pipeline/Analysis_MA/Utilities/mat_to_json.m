% MatToJsonConverter.m
% This script converts a .mat file containing an analysis struct to a JSON file,
% with specific handling for cell arrays containing ConvexHull structures.

% Open a file selection dialog
[fileName, filePath] = uigetfile('*.mat', 'Select the analysis .mat file');

fields_to_take = {
                 'good_frames_animal_identifier',...
                 'annot_reordered',...
                  };
% Check if a file was selected
if fileName == 0
    disp('No file selected. Exiting...');
    return;
end

% Construct the full file path
matFileName = fullfile(filePath, fileName);

% Load the .mat file
data = load(matFileName);

% Assuming the analysis struct is the first variable in the loaded data
structName = fieldnames(data);
analysisStruct = data.(structName{1});

% Function to convert MATLAB data to JSON-compatible format
function jsonCompatible = convertToJsonCompatible(data)
    if iscell(data)
        jsonCompatible = cell(size(data));
        for i = 1:numel(data)
            if isstruct(data{i}) && isfield(data{i}, 'ConvexHull') 
                % Handle ConvexHull structure
                jsonCompatible{i} = convertConvexHull(data{i});
            else
                jsonCompatible{i} = convertToJsonCompatible(data{i});
            end
        end
    elseif isstruct(data)
        jsonCompatible = struct();
        fields = fieldnames(data);
        for i = 1:numel(fields)
            jsonCompatible.(fields{i}) = convertToJsonCompatible(data.(fields{i}));
        end
    elseif ismatrix(data) && ~ischar(data)
        if numel(data) == 1
            jsonCompatible = data;  % Keep scalar values as is
        else
            jsonCompatible = num2cell(data);  % Convert non-scalar matrices to cell arrays
        end
    else
        jsonCompatible = data;
    end

    % Ensure JSON compatibility
    try
        jsonCompatible = ensureJsonCompatible(jsonCompatible);
    catch
        disp('!')
    end

end

% Function to convert ConvexHull structure to JSON-compatible format
function jsonConvexHull = convertConvexHull(convexHull)
    jsonConvexHull = struct();
    if isfield(convexHull, 'convexHull')
        jsonConvexHull.ConvexHull = convexHull.ConvexHull;
    end
    if isfield(convexHull, 'SimplexIndices')
        jsonConvexHull.SimplexIndices = convexHull.SimplexIndices;
    end
    % Add other relevant fields of ConvexHull as needed
    % For example:
    % if isfield(convexHull, 'Volume')
    %     jsonConvexHull.Volume = convexHull.Volume;
    % end
end


% Function to ensure JSON compatibility
function data = ensureJsonCompatible(data)
    if isstruct(data)
        fields = fieldnames(data);
        for i = 1:numel(fields)
            data.(fields{i}) = ensureJsonCompatible(data.(fields{i}));
        end
    elseif iscell(data)
        for i = 1:numel(data)
            data{i} = ensureJsonCompatible(data{i});
        end
    elseif ~(isnumeric(data) || islogical(data) || ischar(data) || isempty(data))
        % Convert any non-standard data to string
        data = char(string(data));
    end
end

% Convert the analysis struct to JSON-compatible format
jsonCompatibleStruct = struct();
fields = fields_to_take;

% Display the number of fields
disp(['Number of fields in the analysis struct: ' num2str(numel(fields))]);

% Iterate through all fields
for i = 1:numel(fields)
    jsonCompatibleStruct.(fields{i}) = convertToJsonCompatible(analysisStruct.(fields{i}));
end

% Convert to JSON with error handling
try
    jsonString = jsonencode(jsonCompatibleStruct, 'PrettyPrint', true);
catch ME
    disp('Error in JSON encoding. Attempting to identify problematic fields...');
    fields = fieldnames(jsonCompatibleStruct);
    for i = 1:numel(fields)
        try
            jsonencode(jsonCompatibleStruct.(fields{i}));
        catch
            disp(['Problematic field: ' fields{i}]);
            % Optionally, you can remove or modify the problematic field
            jsonCompatibleStruct = rmfield(jsonCompatibleStruct, fields{i});
        end
    end
    jsonString = jsonencode(jsonCompatibleStruct, 'PrettyPrint', true);

    rethrow(ME);
end


% Generate output JSON file name
[~, name, ~] = fileparts(fileName);
jsonFileName = fullfile(filePath, [name '_analysis.json']);

% Write JSON to file
fid = fopen(jsonFileName, 'w');
fprintf(fid, '%s', jsonString);
fclose(fid);

disp(['JSON file created: ' jsonFileName]);
new_fields = fieldnames(jsonCompatibleStruct);
disp(['Number of fields processed: ' num2str(numel(new_fields))]);