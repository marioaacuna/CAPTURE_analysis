function M = animate_leg_markers_aligned(mocapstruct, frame_inds, fighand, texthere, markers_to_use)
%matlab_fr = 10;
if nargin < 3
    h = figure(370);
    texthere = '';
else
    h = fighand;
end

frame_last = 0;

% Specify the markers to use
if nargin < 5
    markers_to_use = {'SpineF', 'SpineM', 'KneeR', 'AnkleR', 'HindpawR'};
end
marker_plot = ones(1, numel(markers_to_use));

% Initialize the figure
set(h, 'Color', 'k')

xx = squeeze(mocapstruct.markers_aligned_preproc.(markers_to_use{1})(1, 1));
yy = squeeze(mocapstruct.markers_aligned_preproc.(markers_to_use{1})(1, 2));
zz = squeeze(mocapstruct.markers_aligned_preproc.(markers_to_use{1})(1, 3));
handle_base = line(xx, yy, zz, 'Marker', 'o', 'Color', mocapstruct.markercolor{1}, 'MarkerFaceColor', mocapstruct.markercolor{1}, 'MarkerSize', 2);

ax = gca;
axis(ax, 'manual')
set(gca, 'Color', 'k')
grid on;
set(gca, 'Xcolor', [1 1 1]);
set(gca, 'Ycolor', [1 1 1]);
set(gca, 'Zcolor', [1 1 1]);

zlim([-50 50])
xlim([-30 40])
ylim([-50 50])

set(gca, 'XTickLabels', [], 'YTickLabels', [], 'ZTickLabels', [])
view([-0, 0]);
links_to_use =[5,20,21,22];
for lk = reshape(frame_inds, 1, [])
    cla;

    ind_to_plot = lk;

    % Plot markers that are tracked in the frame
    set(gca, 'Nextplot', 'ReplaceChildren');
    handles_here = cell(1, numel(markers_to_use));
    title(texthere, 'Color', 'w')
    for jj = 1:numel(markers_to_use)
        % Don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(markers_to_use{jj})(ind_to_plot, :), 2))
            if (~sum(mocapstruct.markers_preproc.(markers_to_use{jj})(ind_to_plot, :), 2) == 0)
                xx = squeeze(mocapstruct.markers_aligned_preproc.(markers_to_use{jj})(ind_to_plot, 1));
                yy = squeeze(mocapstruct.markers_aligned_preproc.(markers_to_use{jj})(ind_to_plot, 2));
                zz = squeeze(mocapstruct.markers_aligned_preproc.(markers_to_use{jj})(ind_to_plot, 3));
                handles_here{jj} = line(xx, yy, zz, 'Marker', 'o', 'Color', mocapstruct.markercolor{jj}, 'MarkerFaceColor', mocapstruct.markercolor{jj}, 'MarkerSize', 15);

                hold on
                marker_plot(jj) = 1;
            else
                marker_plot(jj) = 0;
            end
        end
    end

    % Plot the links between markers
    for mml = 1:numel(links_to_use)
        mm = links_to_use(mml);
        % if numel(links_to_use(mm))
            % if (ismember(links_to_use(mm)(1), 1:numel(markers_to_use)) && ismember(links_to_use(mm)(2), 1:numel(markers_to_use)))
                % if (marker_plot(mocapstruct.links{mm}(1)) == 1 && marker_plot(mocapstruct.links{mm}(2)) == 1)
                    xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot, 1)) ...
                          squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot, 1))];
                    yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot, 2)) ...
                          squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot, 2))];
                    zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot, 3)) ...
                          squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot, 3))];
                    line(xx, yy, zz, 'Color', mocapstruct.markercolor{mocapstruct.links{mm}(1)}, 'LineWidth', 3);
                % end
            % end
        % end
    end
    title(lk)
    drawnow
    hold off

    frame_last = lk;

    M(find(frame_inds == lk)) = getframe(gcf);
end
end
