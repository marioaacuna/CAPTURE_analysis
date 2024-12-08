% analysis of movement

% Initialization
logger('Starting movement analysis script', 'INFO');
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
logger('Loading data', 'INFO');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');

% Load analysis structre
load(GC.filename_analysis, 'analysisstruct');

load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;


% Preprocess data
logger('Preprocessing data', 'INFO');
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 15);

% run movement pattern function
logger('Running movement pattern function', 'INFO');
m_f = analyze_movement_patterns(markers_aligned_ds);

% Extract conditions
logger('Extracting conditions', 'INFO');
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

% Run statistics across experimental conditions
logger('Running statistics across experimental conditions', 'INFO');
fields = fieldnames(m_f);
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
          [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};

for f = 1%:length(fields)

    fieldname = fields{f};
    subfields = fieldnames(m_f.(fieldname));
    
    for s = 1:length(subfields)
        subfieldname = subfields{s};
        if startsWith(subfieldname, 'coh'), continue,end
        data = m_f.(fieldname).(subfieldname);

        baseline_data = data(ismember(conditions,'B'),:);
        % Initialize cell array to store data for each condition
        condition_data = cell(length(unique_conditions), 1);
        
        for c = 1:length(unique_conditions)
            condition = unique_conditions{c};
            condition_mask = strcmp(conditions, condition);
            condition_data{c} = data(condition_mask, :) - mean(baseline_data);
        end
        
        % Perform statistical analysis
        all_data = [];
        group_labels = [];
        for c = 1:length(unique_conditions)
            all_data = [all_data; condition_data{c}];
            group_labels = [group_labels; repmat(unique_conditions(c), size(condition_data{c}, 1), 1)];
        end
        
        % Perform Kruskal-Wallis test
        [p_kw, tbl_kw, stats_kw] = kruskalwallis(all_data(:,3), group_labels, 'off');
        
        % Perform ANOVA test
        [p_anova, tbl_anova, stats_anova] = anova1(all_data(:,3), group_labels, 'off');
        
        % Display results
        fprintf('Kruskal-Wallis test for %s - %s, p-value: %.4f\n', fieldname, subfieldname, p_kw);
        fprintf('ANOVA test for %s - %s, p-value: %.4f\n', fieldname, subfieldname, p_anova);
        
        % Create figure for comparison
        figure('Position', [100 100 800 600], 'Color', 'w');
      
        % Calculate means and SEMs for each condition
        means = zeros(1, length(unique_conditions));
        sems = zeros(1, length(unique_conditions));
        for c = 1:length(unique_conditions)
            means(c) = median(condition_data{c}(:, 3));  % z plane
            sems(c) = std(condition_data{c}(:, 3)) / sqrt(size(condition_data{c}, 1));
        end

        % Create bar plot with error bars
        b = bar(means, 'FaceColor', 'flat');
        hold on;
        errorbar(1:length(unique_conditions), means, sems, 'k', 'LineStyle', 'none', 'CapSize', 10);
        
        % Customize plot
        for c = 1:length(unique_conditions)
            b.CData(c,:) = colors{c};
        end
        
        xlabel('Condition');
        ylabel(subfieldname);
        ylim([-0.2 0.5])
        title(sprintf('%s - %s Across Conditions', fieldname, subfieldname), 'Interpreter','none');
        set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
        box off;
        set(gca, 'TickDir', 'out');

        % Add significance marker if test is significant
        if p_kw < 0.05
            plot(1:length(unique_conditions), max(means + sems) * 1.1 * ones(1,length(unique_conditions)), 'k-');
            text(mean(1:length(unique_conditions)), max(means + sems) * 1.15, sprintf('KW p = %.3f', p_kw), ...
                'HorizontalAlignment', 'center');
        end
        if p_anova < 0.05
            plot(1:length(unique_conditions), max(means + sems) * 1.2 * ones(1,length(unique_conditions)), 'k--');
            text(mean(1:length(unique_conditions)), max(means + sems) * 1.25, sprintf('ANOVA p = %.3f', p_anova), ...
                'HorizontalAlignment', 'center');
        end
    end
end

disp('done')
logger('Movement analysis script completed', 'INFO');

    % Negative values indicate guarding (affected side moves less)
    % Positive values indicate compensation (affected side moves more)