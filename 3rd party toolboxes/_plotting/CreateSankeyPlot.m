function CreateSankeyPlot(InputData)
% Written by roman.kuster@alumni.ethz.ch (24.06.2021)
% Function to create a sankey plot / alluvial flow diagram with numerical
% data stored in InputData (e.g. categories 0, 1, and 2). Each row in 
% InputData represents a subject, and each column represents a timepoint.
% The entries in InputData reflect the group a particular subject was
% assigned to at a particular timepoint, and the plot visualises the flow
% of group membership over all timepoints.
%
% - size(InputData,1) corresponds to the number of subjects
% - size(InputData,2) corresponds to the number of timepoints
% - entries represent the category a subject belongs to at a given timepoint
%
% Example Data with 5 subjects, 4 timepoints, and 3 categories
% InputData = [0,0,1,0;2,1,2,2;1,1,0,1;0,1,1,2;1,2,2,2];
% Te highest category (2 in example data) is assigned to the top bar, the
% lowest category to the bottom bar. Can be modified by change the order in
% categories. The bars are coloured from green (top) to red (bottom), and 
% the percentage of subjects belonging to each category is indicated in the
% middle of the bar.
%
% Function calls the modified function sankey_alluvialflow initially
% written by Alexander Carmeli (Alluvial flow diagram, available from: 
% www.mathworks.com/matlabcentral/fileexchange/66746-alluvial-flow-diagram)
% and modified by Ranran Wang (Sankey Diagram, available from:
% www.mathworks.com/matlabcentral/fileexchange/75813-sankey-diagram).
%% Detect categories used:
categories = sort(unique(InputData),'descend');
%% stacked Bars for each category at each timepoint t:
for t = 1:size(InputData,2)
    for c = 1:length(categories)
        Bars{t}(1,c) = sum(InputData(:,t)==categories(c)); % sum categories...
        timePointNames{t} = ['Time ',num2str(t)];
    end
end
%% Change between timepoint t and t+1 for each category:
for t = 1:size(InputData,2)-1
    for c1 = 1:length(categories) % change from 1 category at t...
        for c2 = 1:length(categories) % to all other categories at t+1
            Change{t}(c1,c2) = sum(InputData(:,t) == categories(c1) & InputData(:,t+1) == categories(c2));
        end
    end
end
%% Define figure appearance:
% Bar Colors (from green at the top to red at the bottom)
BarColors(1:length(categories),1) = [0.25: 0.5/(length(categories)-1) :0.75]; % red
BarColors(1:length(categories),2) = [0.75: -0.5/(length(categories)-1) :0.25]; % green
BarColors(1:length(categories),3) = repmat(0.15,length(categories),1); % blue
% Transparency (FaceAlpha) of change-flow:
ChangeTransparancy = 0.2;
% Timepoints on horizontal axis:
X = 0:size(InputData,2)-1;
% Width of the bars:
BarWidth = 30;
%% Create the Figure
figure();
y1_category_points=[];
for t=1:size(InputData,2)-1
    y1_category_points=sankey_alluvialflow(Bars{t}, Bars{t+1}, Change{t}, X(t), X(t+1), y1_category_points,BarColors,ChangeTransparancy,BarWidth);
end
%% Add text (% in each category and timepoint)
ymax = 1.1 * size(InputData,1);
gap = (ymax-size(InputData,1)) / (length(categories)-1);
for t = 1:size(InputData,2)
    for c = 1:length(categories)
        NumbToPrint = num2str(roundn(Bars{t}(c)/size(InputData,1)*100,0));
        text(X(t), Bars{t}(c)/2 + sum(Bars{t}(1:c-1)) + (c-1)*gap, [NumbToPrint,'%'],'HorizontalAlignment','center','Color',[0 0 0])
    end
    text(X(t), Bars{t}(c) + sum(Bars{t}(1:c-1)) + (c+1)*gap, timePointNames{t},'HorizontalAlignment','center','Interpreter','None')
end
return
