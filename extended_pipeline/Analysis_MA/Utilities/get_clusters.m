function clusters = get_clusters(animal_list, animal_ID, cond_inds, analysisstruct, upsamplig_factor, conditions)
    % 1. Find index of animal '328'
    animal_idx = find(strcmp(animal_list, animal_ID));

    % 2. Get the condition indices for animal '328'
    specific_cond_inds = find(ismember(cond_inds,animal_idx));

    % 3. Extract specific frames and outcomes based on the condition indices
    % analysisstruct.frames_with_good_tracking is at least diff 50 frames
    frames = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds);
    outcomes = analysisstruct.annot_reordered_matched{1, 1}(specific_cond_inds); % sampled every 50
    % analysisstruct.frames_tracking_appendages(specific_cond_inds) % 
    % get total number of frames: for this we need to search for the prediction.mat file in the animal folder
    % load the prediction file
    % load(fullfile("D:\test_CAPTURE",animal_ID, 'predictions.mat' ), "predictions");
    this_exp_cond = conditions{animal_idx};
    % if strcmp(this_exp_cond, 'F')
    %     folder_exp_cond = 'PFA';
    % elseif strcmp(this_exp_cond, 'S')
    %     folder_exp_cond = 'saline';
    % end

    disp(['Running Animal ', animal_ID, ' - condition ', this_exp_cond])
    % server_folder = fullfile("H:\DANNCE\6cam_behavior",folder_exp_cond, animal_ID,"DANNCE_ready\DANNCE\predict_results_net_8" );
    % use com to read the number of frames
    % TODO: check if this is correct, becuase for some reasong for an example animal com has 216001 frames,
    % but the predictions has 216000
    % com_filename = fullfile(server_folder, 'com3d_used.mat');
    % COM = load(com_filename, "com");
    n_frames_total = 30*60*100; % TODO, fix later
    % n_frames_total = length(frames);


    % Initialize the clusters array
    clusters = zeros(1, n_frames_total * upsamplig_factor);

    for i = 1:length(frames)
        % frame_val = ceil(frames_this_idx(i) / upsamplig_factor);

        frame_val = frames(i);

        if i == 1
            relative_start = frame_val; % This will be adjusted below
            if animal_idx > 1 && frame_val > analysisstruct.tsnegranularity  % 50
                previous_frame = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds(1) - 1); % diff of 50 frames 
                % relative_start = frame_val - ceil(previous_frame / 3);
                relative_start = frame_val - ceil(previous_frame);
                if relative_start == analysisstruct.tsnegranularity % means it's only one frame
                    relative_start = 1;
                elseif relative_start > analysisstruct.tsnegranularity % meaning at the beginning there are no clusters
                    % relative_start = previous_frame + 50; TODO, still I don't
                    % know how to go on with this
                    keyboard %#ok<KEYBOARDFUN>
                end
            end
            clusters(1:relative_start) = outcomes(i);
            % clusters(relative_start:frame_val) = outcomes(i);

            current_index = relative_start + 1;
        else

            % set the gap: cehck if the previous frame has a gap of 50
            gap = ceil((frames(i) - frames(i-1)));
            if gap > ceil(analysisstruct.tsnegranularity) 
                % keyboard % TODO
                % clusters(current_index:current_index+gap-1) = 0; % fill gaps with zeros
                clusters(current_index:current_index+gap-1) = NaN; % Fill the gaps with NaN. TODO: I still don't know why we have these gaps (check script 1) 
                current_index = current_index + gap;
            else
                clusters(current_index:current_index+gap-1) = outcomes(i);
                % current_index = current_index + gap-1;
                current_index = current_index + gap;

            end
        end
    end
    % Downsample the clusters
    clusters = downsample(clusters, upsamplig_factor);
end