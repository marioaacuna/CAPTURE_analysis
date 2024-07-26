% This file contains a set of general configurations. All variables should
% be loaded into a structure called GC.

function GC = general_configs()
    % Initialize structure
    GC = struct();
    
    % Set version of repository
    GC.version = '0.0.1';
    %% Get path of this file
    current_path = mfilename('fullpath');
    % Remove filename to get root path of the repository
    repository_root_path = regexp(current_path, filesep(), 'split');
    GC.repository_root_path = fullfile(repository_root_path{1:end-1});
    % Temp root folder for outputs
    if ispc
        temp_root = 'D:/CAPTURE';
        project_path = 'C:\Users\acuna\OneDrive - Universitaet Bern\Spontaneous_pain_kinematics\';
    else
        mac_name = 'marioacuna';
        project_path =  fullfile('Users',mac_name,'Library/CloudStorage/OneDrive-UniversitaetBern/Spontaneous_pain_kinematics');
        if ~exist(project_path, 'dir')
            mac_name = 'mario';
            project_path =  fullfile('Users',mac_name,'Library/CloudStorage/OneDrive-UniversitaetBern/Spontaneous_pain_kinematics');
        end

        temp_root = fullfile('/Users',mac_name,'Documents/Temp_analysis/CAPTURE');
    end
    if ~exist(temp_root, 'dir')
        mkdir(temp_root)
    end
    % Preprocessing folder: Where to store the concatenated files, analysis struct, etc.
    preprocessing_rootpath = fullfile(project_path,'data','0_preprocessing');
    if ~exist(preprocessing_rootpath, 'dir')
        mkdir(preprocessing_rootpath)
    end

    % Preprocessing folder: Where to store the concatenated files, analysis struct, etc.
    postprocessing_rootpath = fullfile(project_path,'data','1_postprocessing');
    if ~exist(postprocessing_rootpath, 'dir')
        mkdir(postprocessing_rootpath)
    end

    % Figure folder
    figure_folder = fullfile(project_path, 'figures');
    if ~exist(figure_folder, 'dir')
        mkdir(figure_folder)
    end

    %% Allocate variables
    GC.temp_root                = temp_root; 
    GC.data_rootpath            = data_rootpath;
    GC.preprocessing_rootpath   = preprocessing_rootpath;
    GC.postprocessing_rootpath  = postprocessing_rootpath;
    GC.figure_folder            = figure_folder;

    % Python
    GC.python = struct();
    GC.python.environment_name = 'base';
    [~, msg] = system(sprintf('activate %s && python -c "import sys; print(sys.executable)"', GC.python.environment_name));
    GC.python.interpreter_path = msg(1:end-1);
    GC.python.scripts_path = fullfile(GC.repository_root_path, 'Utilities', 'python');
    
    % R
    GC.R = struct();
    GC.R.scripts_path = fullfile(GC.repository_root_path, 'Code', 'R');

    
    %% PREPROCESSING   
   
    %% PLOTS
   
    % Set plotting options for graphs made in python
    GC.python.font_size_labels  = 16;
    GC.python.scatterplot_small = 7;
    GC.python.scatterplot_large = 10;
    
    
    
    %% TO BE FILLED IN GUI
    GC.experiment_name = '';
    

%% MLint exceptions
%#ok<*CTCH>
