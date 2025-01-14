% Initialize containers for per-cluster analysis
H_animals = fieldnames(data_H);
N_animals = fieldnames(data_N);

% Let's first analyze one animal to understand cluster patterns
example_animal = H_animals{1};
num_clusters = length(data_H.(example_animal).ensemble_activity);
num_neurons = size(data_H.(example_animal).ensemble_activity(5).principal_components, 1);

fprintf('Number of clusters: %d\n', num_clusters);
fprintf('Number of neurons: %d\n', num_neurons);

% Create analysis function for ensemble patterns
function analyze_cluster_ensembles(data_struct, animal_id)
    ensemble_data = data_struct.(animal_id).ensemble_activity;
    
    % Initialize matrix to store ensemble weights per cluster
    num_clusters = length(ensemble_data);
    num_neurons = size(ensemble_data(5).neural_weights, 1);
    cluster_weights = zeros(num_clusters, num_neurons);
    
    % Fill matrix with weights
    for cluster = 1:num_clusters
        if ~isempty(ensemble_data(cluster).neural_weights)
            cluster_weights(cluster, :) = ensemble_data(cluster).neural_weights';
        end
    end
    
    % Visualize ensemble patterns
    figure('Position', [100 100 1200 600]);
    
    % 1. Heatmap of neural weights across clusters
    subplot(1,2,1)
    imagesc(cluster_weights)
    colorbar
    title(['Neural Ensembles Across Clusters - Animal ' animal_id])
    xlabel('Neuron ID')
    ylabel('Cluster ID')
    colormap('jet')
    
    % 2. Find consistently co-active neurons
    correlation_matrix = corr(cluster_weights');
    subplot(1,2,2)
    imagesc(correlation_matrix)
    colorbar
    title('Cluster Similarity Based on Neural Ensembles')
    xlabel('Cluster ID')
    ylabel('Cluster ID')
    colormap('jet')
    
    % Return important metrics
    return_data.cluster_weights = cluster_weights;
    return_data.correlation_matrix = correlation_matrix;
    
    % Find highly consistent ensembles (clusters with similar patterns)
    threshold = 0.7; % Correlation threshold
    [rows, cols] = find(correlation_matrix > threshold & correlation_matrix < 1);
    
    fprintf('\nHighly similar clusters (r > %.2f):\n', threshold);
    for i = 1:length(rows)
        if rows(i) < cols(i)
            fprintf('Clusters %d and %d: r = %.2f\n', rows(i), cols(i), correlation_matrix(rows(i), cols(i)));
        end
    end
    
    % Identify most consistently active neurons
    mean_weights = mean(cluster_weights, 1);
    [sorted_weights, sorted_idx] = sort(mean_weights, 'descend');
    
    fprintf('\nTop 5 most consistently active neurons:\n');
    for i = 1:5
        fprintf('Neuron %d: Mean weight = %.3f\n', sorted_idx(i), sorted_weights(i));
    end
end

% Analyze example animals from each condition
analyze_cluster_ensembles(data_H, H_animals{1});
title('Sham Condition')

analyze_cluster_ensembles(data_N, N_animals{1});
title('Neuropathic Condition')

%% Concatenated analysis of PCA data
% Initialize group_analysis structure
group_analysis = struct();
group_analysis.H = struct();
group_analysis.N = struct();

% Get animal IDs
H_animals = fieldnames(data_H);
N_animals = fieldnames(data_N);

% First, let's understand the structure dimensions
fprintf('Number of Sham animals: %d\n', length(H_animals));
fprintf('Number of Neuropathic animals: %d\n', length(N_animals));

% Let's look at one animal to understand the structure
example_animal = H_animals{1};
fprintf('Example animal structure:\n');
disp(size(data_H.(example_animal).ensemble_activity(5).principal_components));
disp(size(data_H.(example_animal).ensemble_activity(5).neural_weights));
% Initialize arrays
H_components = [];
H_weights = [];
H_clusters = [];
H_pca_scores = [];

% For Sham animals
for h = 1:length(H_animals)
    animal_ID = H_animals{h};
    animal_data = data_H.(animal_ID).ensemble_activity;
    
    fprintf('Processing animal %s\n', animal_ID);
    fprintf('Number of clusters: %d\n', length(animal_data));
    
    % For each pose cluster in this animal
    for cluster = 1:length(animal_data)
        % Check if this cluster has activity (non-empty components and weights)
        if ~isempty(animal_data(cluster).principal_components) && ...
           ~isempty(animal_data(cluster).neural_weights)
            
            [rows_pc, cols_pc] = size(animal_data(cluster).principal_components);
            [rows_w, cols_w] = size(animal_data(cluster).neural_weights);
            
            fprintf('Cluster %d: PC size [%d,%d], Weights size [%d,%d]\n', ...
                    cluster, rows_pc, cols_pc, rows_w, cols_w);
            
            % Only concatenate if we have data
            H_components = [H_components; animal_data(cluster).principal_components];
            H_weights = [H_weights; animal_data(cluster).neural_weights];
            H_pca_scores = [H_pca_scores; animal_data(cluster).pca_score];
            
        else
            fprintf('Cluster %d: No activity\n', cluster);
        end
        
    end
    H_clusters = [H_clusters; data_H.(animal_ID).cluster_vector'];
end

% Do the same for Neuropathic animals
N_components = [];
N_weights = [];
N_clusters = [];
N_pca_scores = [];

for n = 1:length(N_animals)
    animal_ID = N_animals{n};
    animal_data = data_N.(animal_ID).ensemble_activity;
    
    fprintf('Processing animal %s\n', animal_ID);
    fprintf('Number of clusters: %d\n', length(animal_data));
    
    for cluster = 1:length(animal_data)
        if ~isempty(animal_data(cluster).principal_components) && ...
           ~isempty(animal_data(cluster).neural_weights)
            
            N_components = [N_components; animal_data(cluster).principal_components];
            N_weights = [N_weights; animal_data(cluster).neural_weights];
            N_pca_scores = [N_pca_scores; animal_data(cluster).pca_score];
        end
    end
    N_clusters = [N_clusters; data_N.(animal_ID).cluster_vector'];

end

% Final size check
fprintf('\nFinal sizes:\n');
fprintf('Sham Components: %d x %d\n', size(H_components));
fprintf('Sham Weights: %d x %d\n', size(H_weights));
fprintf('Neuropathic Components: %d x %d\n', size(N_components));
fprintf('Neuropathic Weights: %d x %d\n', size(N_weights));

%%
% Store in group_analysis
group_analysis = struct();
group_analysis.H.components = H_components;  % Should be Nx3 where N = total clusters * 71
group_analysis.H.weights = H_weights;        % Should be Nx1
group_analysis.H.clusters = H_clusters;        % Should be Nx1
group_analysis.N.components = N_components;
group_analysis.N.weights = N_weights;
group_analysis.N.clusters = N_clusters;
group_analysis.H.pca_scores = H_pca_scores;
group_analysis.N.pca_scores = N_pca_scores;

% Create visualization
figure('Position', [100 100 1200 800]);


% 1. Neural Weights Distribution by Condition
subplot(2,2,1)
% Calculate mean and SEM
mean_H = mean(group_analysis.H.weights);
mean_N = mean(group_analysis.N.weights);
sem_H = std(group_analysis.H.weights)/sqrt(length(group_analysis.H.weights));
sem_N = std(group_analysis.N.weights)/sqrt(length(group_analysis.N.weights));

% Create bar plot
b = bar([1 2], [mean_H mean_N]);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1]; % Blue for Sham
b.CData(2,:) = [1 0 0]; % Red for Neuropathic

hold on
% Add error bars
errorbar([1 2], [mean_H mean_N], [sem_H sem_N], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize plot
set(gca, 'XTick', 1:2)
set(gca, 'XTickLabel', {'Sham', 'Neuropathic'})
title('Neural Weights by Condition')
ylabel('Mean Weight Value ± SEM')

% Add statistical comparison
[~,p] = ttest2(group_analysis.H.weights, group_analysis.N.weights);
text(1.5, max([mean_H+sem_H, mean_N+sem_N])*1.1, sprintf('p = %.3f', p))

%% 2. Component Space Comparison
subplot(2,2,2)
% Create 3D scatter plot
scatter3(group_analysis.H.pca_scores(:,1), ...
        group_analysis.H.pca_scores(:,2), ...
        group_analysis.H.pca_scores(:,3), ...
        20, 'b.', 'DisplayName', 'Sham')
hold on
scatter3(group_analysis.N.pca_scores(:,1), ...
        group_analysis.N.pca_scores(:,2), ...
        group_analysis.N.pca_scores(:,3), ...
        20, 'r.', 'DisplayName', 'Neuropathic')

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('Neural Activity in PC Space')
legend('Location', 'best')
grid on
view(45, 30)  % Set viewing angle

% Optional: Add rotation capability
rotate3d on
%% 3. Component Correlation Matrix
subplot(2,2,3)
corr_matrix = zeros(3,3);
for i = 1:3
    for j = 1:3
        corr_matrix(i,j) = corr(group_analysis.H.components(:,i), ...
                               group_analysis.N.components(:,j));
    end
end
imagesc(corr_matrix)
colorbar
title('Component Correlation Between Conditions')
xlabel('Neuropathic PCs')
ylabel('Sham PCs')
colormap('jet')
axis square

%% 4. Cumulative Variance by Component
subplot(2,2,4)
% Calculate explained variance for each condition
var_H = var(group_analysis.H.components);
var_N = var(group_analysis.N.components);
cumvar_H = cumsum(var_H)/sum(var_H) * 100;
cumvar_N = cumsum(var_N)/sum(var_N) * 100;

plot(1:3, cumvar_H, 'b-o', 'LineWidth', 2, 'DisplayName', 'Sham')
hold on
plot(1:3, cumvar_N, 'r-o', 'LineWidth', 2, 'DisplayName', 'Neuropathic')
xlabel('Principal Component')
ylabel('Cumulative Variance Explained (%)')
title('Cumulative Variance Explained')
legend('Location', 'southeast')
grid on

% Print summary statistics
fprintf('\nSummary Statistics:\n')
fprintf('Number of samples:\n')
fprintf('  Sham: %d\n', size(group_analysis.H.weights,1))
fprintf('  Neuropathic: %d\n', size(group_analysis.N.weights,1))
fprintf('\nWeight comparison p-value: %.3f\n', p)