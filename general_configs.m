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
    GC.project_path             = project_path;
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
    % file names
    GC.filename_analysis    = fullfile(preprocessing_rootpath,'raw_concat_analysis.mat' );
    GC.filename_predictions = fullfile(preprocessing_rootpath, 'agg_predictions.mat');
    GC.filename_ratception  = fullfile(preprocessing_rootpath, 'ratception_prediction.mat');
    % Vars
    GC.frame_rate           = 100; % Original video recordings FR  
    GC.ratception_name      = 'mario_mouse22';
    GC.linkname             = 'mario_mouse22';
    GC.repfactor            = 3; % Repetition factor for the video; round(300/init_frame_rate)
    GC.tsnegranularity      = 25;
    GC.perplexity           = 200;

    % Clustering
    GC.density_res          = 1001; %resolution of the map
    GC.density_width        = 1;%this depends on the number of frames and how data is clustered, 1 works well.
    GC.expansion_factor     = 1.1; %add a little room to the map after kernel smoothing
    GC.density_threshold    = 1*10^(-5); %remove regions in plots with low density
    %% PLOTS
    
    % Set plotting options for graphs made in python
    GC.python.font_size_labels  = 16;
    GC.python.scatterplot_small = 7;
    GC.python.scatterplot_large = 10;

    
    
    
    %% TO BE FILLED IN GUI
    GC.experiment_name = '';
    

%% MLint exceptions
%#ok<*CTCH>
