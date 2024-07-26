% load mocap

mouse_markers = mocapstruct.markers_aligned_preproc;

% Extract fieldnames from the structure
marker_fields = fieldnames(mouse_markers);

% Initialize an empty matrix
all_markers = [];

% Loop over each field in the structure
for i = 1:numel(marker_fields)
    % Get the data for the current marker
    curr_marker = mouse_markers.(marker_fields{i});
    
    % % Reshape the current marker data to 1D if it's not
    % if size(curr_marker,2) ~= 1
    %     curr_marker = reshape(curr_marker, [], 1);
    % end
    % 
    % Concatenate the current marker data to the all_markers matrix
    all_markers = [all_markers, curr_marker]; 
end

%% standarize
% we need to transpose the matrix
mouse_keypoints = all_markers';
mouse_keypoints_stand = zscore(mouse_keypoints);

%% run pca
[~,score] = pca(mouse_keypoints_stand', 'NumComponents',10);


score = score';
figure, imagesc(score)
%% compute wavelet transform
widths = linspace(10, 100, 15);

% MATLAB's cwt requires scales instead of widths. As a rule of thumb, scale = width / sqrt(2)
scales = widths / sqrt(2);

% First row of the score matrix
score_init = score(1, :);

% Compute the CWT with the 'mexh' (Ricker/Mexican hat wavelet) wavelet
wavelets_init = cwt(score_init, scales, 'mexh');

% Assuming wavelets_init already has the CWT coefficients for the first row

% Loop over rows of score from the second row onwards
for x = 2:10
    % Perform the CWT on the current row of score
    cwtmatr = cwt(score(x, :), scales, 'mexh');
    
    % Concatenate the current CWT coefficients to wavelets_init
    wavelets_init = [wavelets_init; cwtmatr];
end

% Standardize wavelets_init
wavelets_stand = zscore(wavelets_init);

[~,score_wavelet] = pca(wavelets_stand', 'NumComponents',10);
score_wavelet = score_wavelet';


%% calculate the tsne
tsne_em_pose = tsne(score(:,1:100:end)', 'Perplexity',30, 'Verbose',1, 'NumPCAComponents',2, 'NumDimensions',2);
figure, scatter(tsne_em_pose(:,1), tsne_em_pose(:,2), 'bo', 'filled')
title('TSNE pose')

tsne_em_dynamics = tsne(score_wavelet(:,1:100:end)', 'Perplexity',30, 'Verbose',1, 'NumPCAComponents',2, 'NumDimensions',2);
figure, scatter(tsne_em_dynamics(:,1), tsne_em_dynamics(:,2), 'ro', 'filled')
title('TSNE dynamics')

% combo dynamcis
combo_features = [score, score_wavelet];
tsne_em_both = tsne(combo_features(:,1:100:end)', 'Perplexity',30, 'Verbose',1, 'NumPCAComponents',2, 'NumDimensions',2);
figure, scatter(tsne_em_pose(:,1), tsne_em_pose(:,2), 'go', 'filled')
title('TSNE combo')

