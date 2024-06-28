function h = sankey_alluvialflow(Bars1, Bars2, change, x1, x2,last_category_points,BarColors,ChangeTransparancy,BarWidth)
Height = sum(Bars1)*1.1;
Gap = sum(Bars1)*0.1 / (length(Bars1)-1);
axis ij % origin is top left
axis off
    
hold on
    
% These are the top points for each left category, with gaps added.
if isempty(last_category_points)
    y1_category_points = [0 cumsum(Bars1)] + (0:numel(Bars1)) .* Gap;
    y1_category_points(end) = [];
else 
    y1_category_points=last_category_points;
end
% These are the top points for each right category, with gaps added.
y2_category_points = [0 cumsum(Bars2)] + (0:numel(Bars2)) .* Gap;
y2_category_points(end) = [];
h=y2_category_points;
     
% Draw the patches, an entire left category at a time
right_columns_so_far = y2_category_points(1:end); % Start at the beginning of each right category and stack as we go.
patches_per_left_category = size(change, 2);
for k_left = 1:size(change, 1) % for each row
    
    % Calculate the coordinates for all the patches split by the
    % Split the left category
    left_patch_points = [0 cumsum(change(k_left, :))] + y1_category_points(k_left);
    patch_top_lefts = left_patch_points(1:end-1);
    patch_bottom_lefts = left_patch_points(2:end);
    
    % Compute and stack up slice of each right category
    patch_top_rights = right_columns_so_far;
    patch_bottom_rights = patch_top_rights + change(k_left, :);
    right_columns_so_far = patch_bottom_rights;
    
    % Plot the patches
    
    % X coordinates of patch corners
    [bottom_curves_x, bottom_curves_y] = get_curves(x1+0.1, patch_bottom_lefts, x2-0.1, patch_bottom_rights);
    [top_curves_x,    top_curves_y]    = get_curves(x2-0.1, patch_top_rights,   x1+0.1, patch_top_lefts);
    X = [ ...
        repmat([x1; x1], 1, patches_per_left_category); % Top left, bottom left
        bottom_curves_x;
        repmat([x2; x2], 1, patches_per_left_category); % Bottom right, top right
        top_curves_x
        ];
    
    % Y coordinates of patch corners
    Y = [ ...
        patch_top_lefts;
        patch_bottom_lefts;
        bottom_curves_y;
        patch_bottom_rights;
        patch_top_rights;
        top_curves_y
        ];
    
    patch('XData', X, 'YData', Y, 'FaceColor', BarColors(k_left,:), 'FaceAlpha', ChangeTransparancy, 'EdgeColor', 'none');
end % for each row
% plot left category bars
for i=1:numel(y1_category_points)
    y1=[y1_category_points; (y1_category_points + Bars1)];
    plot(ones(2, 1)*x1, y1(:,i), 'Color', BarColors(i,:),'LineWidth',BarWidth);
end
hold on
% plot right category bars
for i=1:numel(y2_category_points)
    y2=[y2_category_points; (y2_category_points + Bars2)];
    plot(ones(2, 1)*x2, y2(:,i), 'Color', BarColors(i,:),'LineWidth',BarWidth);
end
    
end % alluvialflow
function [x, y] = get_curves(x1, y1, x2, y2)
% x1, x2: scalar x coordinates of line start, end
% y1, y2: vectors of y coordinates of line start/ends
    Npoints = 15;
    t = linspace(0, pi, Npoints);
    c = (1-cos(t))./2; % Normalized curve
    
    Ncurves = numel(y1);
    y = repmat(y1, Npoints, 1) + repmat(y2 - y1, Npoints,1) .* repmat(c', 1, Ncurves);
    x = repmat(linspace(x1, x2, Npoints)', 1, Ncurves);
end  % get_curve