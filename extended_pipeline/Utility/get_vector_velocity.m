% function [vectorout,vectorstdout,vectorout_accel,vectorstdout_accel]= get_vector_velocity(vectorin,params)
% 
% if params.medfiltorder>1
% tracesmooth = medfilt1(vectorin,params.medfiltorder);
% else
% tracesmooth = vectorin;
% end
% %gfilter = fspecial('gaussian',[50 1], params.gaussorder);
% %tracesmoothed = convn(tracesmooth,gfilter,'same');
% tracesmoothed = tracesmooth;
% endptcond = 'shrink';
% vectorout_vel = cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
% 
% 
% %% try interpolating  
% 
% 
%     vectorout_vel =  interpolate_ends( vectorout_vel,  3*params.difforder_movav);
% %vectorout_vel = cat(1,zeros(params.difforder_movav*1,1),vectorout_vel,zeros(params.difforder_movav*1,1));
% vectorout = movmean(vectorout_vel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
% % vectorout = vectorout((params.difforder_movav+1):(end-params.difforder_movav));
% % 
% % padsize = (numel(vectorin)-numel(vectorout));
% % vectorout = cat(1,vectorout,zeros(padsize,1));
% 
%     vectorout_vel_big =  interpolate_ends( vectorout_vel,  3*params.difforder_movav);
% vectorstdout = movstd(vectorout_vel_big,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
% % vectorstdout = vectorstdout((params.difforder_movav+1):(end-params.difforder_movav));
% 
%   %vectorstdout = cat(1,vectorstdout,zeros(padsize,1));
% 
% a_diff_order = floor(params.difforder);
% vectorout_vel_a = cat(1,zeros((a_diff_order)-2,1),tracesmoothed(a_diff_order:end,1) - tracesmoothed(1:end-a_diff_order+1),1);
% 
% vectorout_accel = cat(1,zeros((a_diff_order)-2,1),vectorout_vel_a(a_diff_order:end,1) - vectorout_vel_a(1:end-a_diff_order+1),1);
% %vectorout_accel = cat(1,zeros(params.difforder_movav*2,1),vectorout_accel,zeros(params.difforder_movav*2,1));
%     vectorout_accel =  interpolate_ends( vectorout_accel,  3*params.difforder_movav);
% 
% vectorout_accel = movmean(vectorout_accel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
%  %vectorout_accel = vectorout_accel((2*params.difforder_movav+1):(end-2*params.difforder_movav));
% 
% vectorstdout_accel = movstd(vectorout_accel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
%  %vectorstdout_accel = vectorstdout_accel((2*params.difforder_movav+1):(end-2*params.difforder_movav));
% 
% %padsize = (numel(vectorin)-numel(vectorout_accel));
% 
% %vectorout_accel = cat(1,vectorout_accel,zeros(padsize,1));
% 
% %padsize = (numel(vectorin)-numel(vectorstdout_accel));
% %vectorstdout_accel = cat(1,vectorstdout_accel,zeros(padsize,1));
% 
% end

function [vectorout,vectorstdout,vectorout_accel,vectorstdout_accel] = get_vector_velocity(vectorin,params)

% Input validation
n_frames = length(vectorin);

% Smooth input if needed
if params.medfiltorder > 1
    tracesmooth = medfilt1(vectorin, params.medfiltorder);
else
    tracesmooth = vectorin;
end

% Add Gaussian smoothing for better results
if params.gaussorder > 0
    gfilter = fspecial('gaussian', [params.gaussorder*6 1], params.gaussorder);
    tracesmooth = conv(tracesmooth, gfilter, 'same');
end

% Compute velocity using finite differences
diff_order = params.difforder;
vectorout_vel = zeros(n_frames, 1);
vectorout_vel(diff_order:end) = (tracesmooth(diff_order:end) - tracesmooth(1:end-diff_order+1)) / diff_order;

% Define window parameters
half_window = floor(params.difforder_movav/2);

% For large windows, use mirror padding instead of interpolation
if params.difforder_movav >= 100
    % Mirror pad the signal to handle edges better
    pad_size = half_window;
    vectorout_vel_padded = [flipud(vectorout_vel(2:pad_size+1)); 
                           vectorout_vel; 
                           flipud(vectorout_vel(end-pad_size:end-1))];
    
    % Apply moving average on padded signal
    vectorout_padded = movmean(vectorout_vel_padded, params.difforder_movav, 'omitnan');
    vectorstdout_padded = movstd(vectorout_vel_padded, params.difforder_movav, 'omitnan');
    
    % Extract the original size
    vectorout = vectorout_padded(pad_size+1:end-pad_size);
    vectorstdout = vectorstdout_padded(pad_size+1:end-pad_size);
    
    % Only NaN the very edges where we don't have enough data
    edge_nan = min(10, floor(diff_order/2)); % Much smaller edge region
    vectorout(1:edge_nan) = NaN;
    vectorstdout(1:edge_nan) = NaN;
    
else
    % For small windows, use the original approach but with better smoothing
    vectorout = movmean(vectorout_vel, params.difforder_movav, 'omitnan', 'Endpoints', 'shrink');
    vectorstdout = movstd(vectorout_vel, params.difforder_movav, 'omitnan', 'Endpoints', 'shrink');
end

% Compute acceleration with same approach
a_diff_order = floor(params.difforder);
vectorout_vel_a = zeros(n_frames, 1);
if a_diff_order < n_frames
    vectorout_vel_a(a_diff_order:end) = (tracesmooth(a_diff_order:end) - tracesmooth(1:end-a_diff_order+1)) / a_diff_order;
end

vectorout_accel_raw = zeros(n_frames, 1);
if 2*a_diff_order < n_frames
    vectorout_accel_raw(a_diff_order:end) = (vectorout_vel_a(a_diff_order:end) - vectorout_vel_a(1:end-a_diff_order+1)) / a_diff_order;
end

% Process acceleration with same padding approach
if params.difforder_movav >= 100
    % Mirror pad
    pad_size = half_window;
    vectorout_accel_padded = [flipud(vectorout_accel_raw(2:pad_size+1)); 
                             vectorout_accel_raw; 
                             flipud(vectorout_accel_raw(end-pad_size:end-1))];
    
    vectorout_accel_padded_smooth = movmean(vectorout_accel_padded, params.difforder_movav, 'omitnan');
    vectorstdout_accel_padded = movstd(vectorout_accel_padded, params.difforder_movav, 'omitnan');
    
    vectorout_accel = vectorout_accel_padded_smooth(pad_size+1:end-pad_size);
    vectorstdout_accel = vectorstdout_accel_padded(pad_size+1:end-pad_size);
    
    % Minimal edge NaN
    edge_nan = min(10, floor(a_diff_order/2));
    vectorout_accel(1:edge_nan) = NaN;
    
else
    vectorout_accel = movmean(vectorout_accel_raw, params.difforder_movav, 'omitnan', 'Endpoints', 'shrink');
    vectorstdout_accel = movstd(vectorout_accel_raw, params.difforder_movav, 'omitnan', 'Endpoints', 'shrink');
end

% Apply additional smoothing filter for very smooth results
if params.difforder_movav >= 100
    % Extra smoothing pass for large windows
    smooth_window = floor(params.difforder_movav/10);
    vectorout = movmean(vectorout, smooth_window, 'omitnan');
    vectorout_accel = movmean(vectorout_accel, smooth_window, 'omitnan');
end

% Final outlier check only for truly extreme values
threshold_factor = 100; % Very high threshold to only catch real artifacts
for output = {vectorout, vectorstdout, vectorout_accel, vectorstdout_accel}
    data = output{1};
    valid_data = data(~isnan(data));
    if ~isempty(valid_data)
        data_median = median(valid_data);
        data_mad = mad(valid_data, 1);
        if data_mad > 0
            extreme_idx = abs(data - data_median) > threshold_factor * data_mad;
            data(extreme_idx) = NaN;
        end
    end
end

end