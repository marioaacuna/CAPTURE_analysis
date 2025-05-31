% Peamble:

%% init
clear 
clc
close all
% location for baseline

bsl_root ='D:\CAPTURE\output_videos\ID_1386\baseline';

% location for formalin
formalin_root = 'D:\CAPTURE\output_videos\ID_1386\Formalin_injection';

% So far we only have bsl
videos = {'beh.mp4', 'FOV.mp4', 'pred.mp4', 'clusters.mp4', 'raw_miniscope.mp4'};
% 'raw_miniscope.mp4' is optional, so we check if it exists
if ~exist(fullfile(bsl_root, 'raw_miniscope.mp4'), 'file')
    videos = videos(1:end-1);  % Remove the last element if it doesn't exist
end

%% Create 2x2 video layout with black background
% Output video settings
output_filename = fullfile(bsl_root, 'combined_2x2_video.mp4');
output_fps = 30;

% Read video objects
video_readers = cell(length(videos), 1);
video_info = struct('Width', [], 'Height', [], 'NumFrames', []);

for i = 1:length(videos)
    video_path = fullfile(bsl_root, videos{i});
    if exist(video_path, 'file')
        video_readers{i} = VideoReader(video_path);
        video_info(i).Width = video_readers{i}.Width;
        video_info(i).Height = video_readers{i}.Height;
        video_info(i).NumFrames = video_readers{i}.NumFrames;
        fprintf('Loaded video %d: %s (%dx%d, %d frames)\n', i, videos{i}, ...
            video_info(i).Width, video_info(i).Height, video_info(i).NumFrames);
    else
        fprintf('Warning: Video file not found: %s\n', video_path);
        video_readers{i} = [];
    end
end

% Remove empty video readers
valid_videos = ~cellfun(@isempty, video_readers);
video_readers = video_readers(valid_videos);
videos = videos(valid_videos);
video_info = video_info(valid_videos);

% Ensure we have at least 4 videos for 2x2 layout, pad with black if needed
num_videos = length(video_readers);
if num_videos < 4
    fprintf('Warning: Only %d videos found. Padding with black frames.\n', num_videos);
end

% Determine output dimensions
% Find the maximum dimensions to ensure all videos fit
max_width = max([video_info.Width]);
max_height = max([video_info.Height]);

% Calculate 2x2 grid dimensions
grid_width = max_width * 2;
grid_height = max_height * 2;

% Determine the minimum number of frames across all videos
min_frames = inf;
for i = 1:length(video_readers)
    if ~isempty(video_readers{i})
        min_frames = min(min_frames, video_readers{i}.NumFrames);
    end
end

if isinf(min_frames)
    error('No valid video files found.');
end

% Create output video writer
output_video = VideoWriter(output_filename, 'MPEG-4');
output_video.FrameRate = output_fps;
open(output_video);

fprintf('Creating 2x2 video layout (%dx%d) with %d frames...\n', grid_width, grid_height, min_frames);

% Process each frame
for frame_idx = 1:min_frames
    % Create black background frame
    combined_frame = zeros(grid_height, grid_width, 3, 'uint8');
    
    % Position videos in 2x2 grid
    positions = [
        1, 1;                           % Top-left
        1, max_width + 1;              % Top-right  
        max_height + 1, 1;             % Bottom-left
        max_height + 1, max_width + 1  % Bottom-right
    ];
    
    for vid_idx = 1:min(4, length(video_readers))
        if ~isempty(video_readers{vid_idx}) && hasFrame(video_readers{vid_idx})
            % Read frame from current video
            current_frame = readFrame(video_readers{vid_idx});
            
            % Resize frame if necessary to fit in the allocated space
            if size(current_frame, 1) ~= max_height || size(current_frame, 2) ~= max_width
                current_frame = imresize(current_frame, [max_height, max_width]);
            end
            
            % Place frame in the appropriate position
            row_start = positions(vid_idx, 1);
            row_end = row_start + max_height - 1;
            col_start = positions(vid_idx, 2);
            col_end = col_start + max_width - 1;
            
            combined_frame(row_start:row_end, col_start:col_end, :) = current_frame;
        end
    end
    
    % Write combined frame to output video
    writeVideo(output_video, combined_frame);
    
    % Display progress
    if mod(frame_idx, 30) == 0 || frame_idx == min_frames
        fprintf('Progress: %d/%d frames (%.1f%%)\n', frame_idx, min_frames, ...
            100 * frame_idx / min_frames);
    end
end

% Close output video
close(output_video);

% Close all video readers
for i = 1:length(video_readers)
    if ~isempty(video_readers{i})
        delete(video_readers{i});
    end
end

fprintf('2x2 video saved to: %s\n', output_filename);

%% Optional: Display final frame preview
if exist('combined_frame', 'var')
    figure('Name', '2x2 Video Layout Preview', 'Color', 'k');
    imshow(combined_frame);
    title('Final Frame of 2x2 Video Layout', 'Color', 'white');
end
