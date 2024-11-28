function metadata = readExperimentMetadata()
    % GC = general_configs();
    % Read JSON file using MATLAB's built-in json reader
    
    json_path = fullfile(pwd, 'extended_pipeline', 'configs', 'metadata.json');
    metadata = jsondecode(fileread(json_path));
end