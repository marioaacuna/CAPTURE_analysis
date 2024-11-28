% FILEPATH: peak_freq_auc.m
animal_list = {'326', '327', '328', '330', '332', '334', '335', '336'};
conditions = {'F', 'F', 'S', 'S', 'S', 'F', 'F', 'F'};
dfof_rootpath = 'V:\Ca_imaging_pain\4_fluorescence_traces';
dfof_suffix = '_raw_deltaF_over_F';
dfof_prefix = 'JH_';

% example filename: JH_326_F_raw_deltaF_over_F.mat
% example filepath: V:\Ca_imaging_pain\4_fluorescence_traces\JH_326_F_raw_deltaF_over_F.mat

% Define the parameters
threshold = 0.5;
min_peak_distance = 10;
min_peak_prominence = 0.1;

% Loop through animals
data = cell(length(animal_list), 1);
for animal_idx = 1:length(animal_list)
    animal = animal_list{animal_idx};
    condition = conditions{animal_idx};
    
    % Load the data
    filename = [dfof_prefix animal dfof_suffix '.mat'];
    filepath = fullfile(dfof_rootpath, filename);
    load(filepath); % Variable called dFF, cell array of cell arrays. We take the second column, which contains ROis x frames
    dfof = dFF(:,2); % session 2

    % initialize the results
    roi_results = cell(size(dfof, 1), 1);
    % Loop through ROIs
    for roi_idx = 1:size(dfof, 1)
        roi_data = dfof{roi_idx, :};
        
        % Calculate the peaks
        [pks, locs] = findpeaks(roi_data, 'MinPeakHeight', threshold, ...
            'MinPeakDistance', min_peak_distance, ...
            'MinPeakProminence', min_peak_prominence);
        
        % Calculate the frequency
        if ~isempty(locs)
            freq = 1 / mean(diff(locs));
        else
            freq = NaN;
        end
        
        % Calculate the AUC
        auc = trapz(roi_data);

        % Store the results for this ROI
        roi_results{roi_idx} = struct('peaks', pks, 'freq', freq, 'auc', auc);

    end
    % Store the data
    data{animal_idx} = roi_results;
end

%% Do stats for condition 
% Initialize results
results = struct('AUC', [], 'Frequency', [], 'PeakHeight', []);

% Separate data based on conditions
data_S = data(strcmp(conditions, 'S'));
data_F = data(strcmp(conditions, 'F'));

% Helper function to extract specific scalar feature from data
extract_scalar_feature = @(data, feature) cellfun(@(roi_list) cellfun(@(roi) roi.(feature), roi_list), data, 'UniformOutput', false);

% Extract scalar features (AUC and frequency)
auc_S_cell = extract_scalar_feature(data_S, 'auc');
auc_F_cell = extract_scalar_feature(data_F, 'auc');

freq_S_cell = extract_scalar_feature(data_S, 'freq');
freq_F_cell = extract_scalar_feature(data_F, 'freq');

% For peaks
extract_peak_feature = @(data) cellfun(@(roi_list) cellfun(@(roi) mean(roi.peaks), roi_list), data, 'UniformOutput', false);

peak_height_S_cell = extract_peak_feature(data_S);
peak_height_F_cell = extract_peak_feature(data_F);

% Convert cells of arrays into single vectors
vec = @(data_cell) cell2mat(cellfun(@(x) x(:), data_cell, 'UniformOutput', false));

auc_S = vec(auc_S_cell);
auc_F = vec(auc_F_cell);

freq_S = vec(freq_S_cell);
freq_F = vec(freq_F_cell);

peak_height_S = vec(peak_height_S_cell);
peak_height_F = vec(peak_height_F_cell);

% T-tests
[~, results.AUC.pval] = ttest2(auc_S, auc_F);
[~, results.Frequency.pval] = ttest2(freq_S, freq_F);
[~, results.PeakHeight.pval] = ttest2(peak_height_S, peak_height_F);

% Convert results to table
results_table = struct2table(results);

% Display the table
disp(results_table);

% TODO: Add plotting and saving of the results
%% Plotting


% Plotting
figure;

% AUC plot
subplot(3, 1, 1);
plot_data_with_bar(auc_S, auc_F, 'AUC', results.AUC.pval);

% Frequency plot
subplot(3, 1, 2);
plot_data_with_bar(freq_S, freq_F, 'Frequency', results.Frequency.pval);

% Peak Height plot
subplot(3, 1, 3);
plot_data_with_bar(peak_height_S, peak_height_F, 'Peak Height', results.PeakHeight.pval);




%% V2
% Plotting
figure;

% AUC plot
subplot(3, 1, 1);
plot_data_with_jitter(auc_S, auc_F, 'AUC');

% Frequency plot
subplot(3, 1, 2);
plot_data_with_jitter(freq_S, freq_F, 'Frequency');

% Peak Height plot
subplot(3, 1, 3);
plot_data_with_jitter(peak_height_S, peak_height_F, 'Peak Height');

%%
% Define a function to plot data with bar and errorbar
% Define a function to plot data with bar and errorbar
function plot_data_with_bar(data_S, data_F, title_str, p)
    % Number of data points
    num_S = length(data_S);
    num_F = length(data_F);

    % Bar plot with mean values
    mean_S = mean(data_S, 'omitnan');
    mean_F = mean(data_F, 'omitnan');
    b1 = bar(1, mean_S, 'FaceAlpha', 0.7);
    b1.FaceColor = 'blue';
    hold on;
    
    b2 = bar(2, mean_F, 'FaceAlpha', 0.7);
    b2.FaceColor = 'red';

    % Error bars with SEM
    SEM_S = std(data_S, 'omitnan') / sqrt(num_S);
    SEM_F = std(data_F, 'omitnan') / sqrt(num_F);
    errorbar([1, 2], [mean_S, mean_F], [SEM_S, SEM_F], 'k.', 'LineWidth', 1, 'Color', 'white');
    
    % Formatting
    xlim([0 3]);
    xticks([1, 2]);
    xticklabels({'S', 'F'});
    ylabel(title_str, 'Color', 'white');
    title([title_str, ' - p =', num2str(p)]    , 'Color', 'white');
    box off

    % Change axis properties
    ax = gca;
    ax.XColor = 'white';
    ax.YColor = 'white';
    ax.Color  = 'black';
    ax.GridColor = 'white';
    ax.MinorGridColor = 'white';

    % Set figure background color
    fig = gcf;
    fig.Color = 'black';

    hold off;
end





% Define a function to plot data with jitter, bar, and errorbar
function plot_data_with_jitter(data_S, data_F, title_str)
    % Number of data points
    num_S = length(data_S);
    num_F = length(data_F);
    
    % Jittered x values for scatter plot
    jitter_amount = 0.2;
    x_S = 1 + (jitter_amount * randn(num_S, 1));
    x_F = 2 + (jitter_amount * randn(num_F, 1));

    % Scatter plot
    scatter(x_S, data_S, 'r', 'jitter', 'on', 'jitterAmount', jitter_amount);
    hold on;
    scatter(x_F, data_F, 'b', 'jitter', 'on', 'jitterAmount', jitter_amount);

    % Bar plot with mean values
    mean_S = mean(data_S, 'omitnan');
    mean_F = mean(data_F, 'omitnan');
    bar([1, 2], [mean_S, mean_F], 'FaceAlpha', 0.3);

    % Error bars with SEM
    SEM_S = std(data_S, 'omitnan') / sqrt(num_S);
    SEM_F = std(data_F, 'omitnan') / sqrt(num_F);
    errorbar([1, 2], [mean_S, mean_F], [SEM_S, SEM_F], 'k.', 'LineWidth', 2);
    
    % Formatting
    xlim([0.5 2.5]);
    xticks([1, 2]);
    xticklabels({'S', 'F'});
    ylabel(title_str);
    title(title_str);
    hold off;
end
