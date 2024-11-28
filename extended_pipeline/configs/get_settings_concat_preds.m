function settings = get_settings_concat_preds()
    % Check if previous settings exist
    settings_file = fullfile(pwd, 'extended_pipeline', 'configs', 'prediction_concat_settings.mat');
    
    if exist(settings_file, 'file')
        prev_settings = load(settings_file);
        use_previous = input('Use previous prediction concatenation settings? (y/n): ', 's');
        if strcmpi(use_previous, 'y')
            settings = prev_settings.settings;
            return;
        end
    end
    
    % Default settings
    defaults = struct(...
        'overwrite_pred_concat', 1, ...
        'overwrite_ratception', 0, ...
        'overwrite_MLmatobjfile', 0, ...
        'overwrite_coefficient', 0, ...
        'overwrite_zvals', 0, ...
        'do_extra_features', 0, ...
        'plot_poses', 1);
    
    % Create settings menu
    fprintf('\nPrediction Concatenation Settings:\n');
    fprintf('Enter new values or press Enter to keep defaults\n');
    
    settings = defaults;
    fields = fieldnames(defaults);
    
    for i = 1:length(fields)
        field = fields{i};
        default_val = defaults.(field);
        
        % Make the prompt more user-friendly
        prompt = strrep(field, '_', ' ');
        prompt = [upper(prompt(1)) lower(prompt(2:end))];
        
        response = input(sprintf('%s? (current: %d) [Enter to keep]: ', prompt, default_val), 's');
        
        if ~isempty(response)
            settings.(field) = str2double(response);
        end
    end
    
    % Ask to save settings for future use
    save_settings = input('Save these prediction concatenation settings for future use? (y/n): ', 's');
    if strcmpi(save_settings, 'y')
        if ~exist(fileparts(settings_file), 'dir')
            mkdir(fileparts(settings_file));
        end
        save(settings_file, 'settings');
        fprintf('Settings saved.\n');
    end
end