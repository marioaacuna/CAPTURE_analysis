% Figure animal density per cluster
% TODO: run it independenlty, not it runs after script__02_02


sample_to_take = 'F'; % 'F', 'S'

if strcmp(sample_to_take, 'F')
    data = cluster_proportions_F;
    clusters_of_interest = predominantF;
    color = 'r';
else
    data = cluster_proportions_S;
    clusters_of_interest = predominantS;
    color = 'b';
end
% plot

% Assuming clusterProportions_F and predominantF are already defined

% Extract only the predominant clusters
predominantClusters = data(:, clusters_of_interest);

% Calculate the number of animals that have each cluster
animalCount = sum(predominantClusters > 0, 1);

% Calculate the proportion of animals that have each cluster
totalAnimals = size(predominantClusters, 1);
animalProportion = animalCount / totalAnimals;

% Create a bar plot
figure('Position',[10 10 1200 500]);
bar(animalProportion, 'FaceColor',color, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Customize the plot
xlabel('Predominant Cluster Index');
ylabel('Proportion of Animals');
title('Proportion of Animals Exhibiting Predominant Clusters in Formalin Condition');
ylim([0 1]);  % Set y-axis limits from 0 to 1
grid off;

% Add text labels on top of each bar
% for i = 1:length(animalProportion)
%     text(i, animalProportion(i), sprintf('%.2f', animalProportion(i)), ...
%          'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% end

% Adjust x-axis ticks and labels
xticks(1:length(clusters_of_interest));
xticklabels(clusters_of_interest);
xtickangle(45);  % Rotate x-axis labels for better readability if needed

% Add a text box with additional information
infoText = sprintf('Total Animals: %d\nTotal Predominant Clusters: %d', ...
                   totalAnimals, length(clusters_of_interest));
annotation('textbox', [0.15, 0.8, 0.2, 0.1], 'String', infoText, 'EdgeColor', 'none');
box off
set(gca, 'TickDir', 'out')

%% Pie chart
% Assuming we have the following variables:
% predominantF - indices of predominant clusters for Formalin condition
% predominantS - indices of predominant clusters for Saline condition
% totalClusters - total number of clusters
totalClusters = length(data);
% Calculate the number of clusters in each category
numPredominantF = length(predominantF);
numPredominantS = length(predominantS);
numPredominantBoth = length(intersect(predominantF, predominantS));
numOtherClusters = totalClusters - (numPredominantF + numPredominantS - numPredominantBoth);

% Prepare data for the pie chart
data = [numPredominantF - numPredominantBoth, ...
        numPredominantS - numPredominantBoth, ...
        numPredominantBoth, ...
        numOtherClusters];

% Create labels
labels = {'Formalin Only', 'Saline Only', 'Both', 'Other Clusters'};

% Create custom colors
colors = [0.8 0.2 0.2;  % Red for Formalin
          0 0.4 0.8;    % Blue for Saline
          0.5 0.5 0.5;  % Gray for Both
          0.9 0.9 0.9]; % Light gray for Others

% Create the pie chart
figure('Position', [100, 100, 800, 600]);
h = pie(data);

% Customize the appearance
colormap(colors);
legend(labels, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Distribution of Predominant Clusters in Formalin and Saline Conditions');

% Add percentage labels to the pie slices
pText = findobj(h, 'Type', 'text');
percentValues = get(pText, 'String');
combinedLabels = strcat(labels, {' '}, percentValues');
pText = findobj(h, 'Type', 'text');
set(pText, {'String'}, combinedLabels');

% Add text annotation for exact numbers
txt = sprintf(['Formalin Only: %d\n' ...
               'Saline Only: %d\n' ...
               'Both: %d\n' ...
               'Other Clusters: %d\n' ...
               'Total Clusters: %d'], ...
    data(1), data(2), data(3), data(4), totalClusters);
annotation('textbox', [0.15, 0.75, 0.2, 0.2], 'String', txt, 'EdgeColor', 'none');

% for some reason the pie chart cannot be copied to illustrator
% so, we export it firts

% Export as PDF
filename = fullfile(GC.figure_folder,'predominant_clusters_distribution_all_frames.pdf');

% Method 1: For MATLAB R2020a and later
if exist('exportgraphics', 'file')
    exportgraphics(gcf, filename, 'ContentType', 'vector');
else
    % Method 2: For earlier versions of MATLAB
    print(gcf, '-dpdf', '-vector', filename);
end

fprintf('Figure exported as %s\n', filename);