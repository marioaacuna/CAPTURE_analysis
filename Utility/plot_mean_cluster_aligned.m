function M = plot_mean_cluster_aligned(mocapstruct,frame_inds,texthere)
%matlab_fr = 10;
if nargin<2
    % h=figure(370);
    texthere = '';
else
    % h=fighand;
end

frame_last = 0;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
% set(h,'Color','k')

% xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
% yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
% zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
% handle_base = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);


% ax = gca;
% axis(ax,'manual')
set(gca,'Color','k')
% grid on;
% set(gca,'Xcolor',[1 1 1 ]);
% set(gca,'Ycolor',[1 1 1]);
% set(gca,'Zcolor',[1 1 1]);
% 
zlim([-20 50])
xlim([-40 50])
ylim([-50 50])


%     zlim([-210 270])
%     xlim([-240 240])
%     ylim([-240 240])
%

% set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
view([-22, 12]);
% fprintf('n frames %i \n',numel(frame_inds));
positions = cell(1, numel(frame_inds));
for lk = reshape(frame_inds,1,[])%1:10:10000
   
    % cla;


    %mocapstruct.links{20} = [];
    %       mocapstruct.links{22} = [];

    ind_to_plot = lk;

    %% Plot markers that are tracked in the frame
    % set(gca,'Nextplot','ReplaceChildren');
    % handles_here = cell(1,numel(mocapstruct.markernames));
    % title(texthere, 'Color','w')
    % Initialize a 3D matrix to store the positions
    positions{lk} = zeros(numel(mocapstruct.markernames), 3);

    % Loop over each marker
    for jj = 1:numel(mocapstruct.markernames)
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
            if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
                xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                % handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',9);

                % Store the positions in the 3D matrix
                positions{lk}(jj, :) = [xx, yy, zz];

                % hold on
                % marker_plot(jj) = 1;
            else
                % marker_plot(jj) = 0;
            end
        end
    end
end

% Initialize a matrix to store the average positions
average_positions = zeros(numel(mocapstruct.markernames), 3);

% Get a logical index of non-empty cells in positions
nonEmptyFrames = ~cellfun(@isempty, positions);
colors= mocapstruct.markercolor;
% Compute the average position for each marker across all frames
for jj = 1:numel(mocapstruct.markernames)
    % Collect the positions for this marker across all non-empty frames
    marker_positions = cell2mat(cellfun(@(pos) pos(jj, :), positions(nonEmptyFrames), 'UniformOutput', false)');
    % average_positions(jj, :) = median(marker_positions, 1);  % average along rows (i.e., over frames)
    average_positions(jj, :) = (marker_positions(1,:,:));  % average along rows (i.e., over frames)

end
hold on

% plot
for ijj = 1:length(average_positions)
    scatter3(average_positions(ijj,1), average_positions(ijj,2), average_positions(ijj,3),3, 'MarkerFaceColor', colors{ijj}, 'MarkerEdgeColor','none')
end
%%

%% plot the links between markers
for mm = 1:numel(mocapstruct.links)

    xx = [squeeze(average_positions(mocapstruct.links{mm}(1),1))...
        squeeze(average_positions(mocapstruct.links{mm}(2),1))];
    yy = [squeeze(average_positions(mocapstruct.links{mm}(1),2))...
        squeeze(average_positions(mocapstruct.links{mm}(2),2))];

    zz = [squeeze(average_positions(mocapstruct.links{mm}(1),3))...
        squeeze(average_positions(mocapstruct.links{mm}(2),3))];


    % xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
    %     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
    % yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
    %     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
    % zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
    %     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
    line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',3);
end
hold off

% drawnow
% hold off
% 
% frame_last = lk;
% 
% M(find(frame_inds == lk)) =  getframe(gcf);
% 
% %clf

%set(gca,'Nextplot','add');

end