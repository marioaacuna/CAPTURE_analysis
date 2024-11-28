function frames = timeToFrames(timeString, FR)
    % Convert time string to frames
    % Input: timeString in format 'MMSSXX' where XX is milliseconds/10
    % Output: frame number at 100 fps, rounded to nearest integer
    % Example usage:
    % timeString = '245668';
    % frameNumber = timeToFrames(timeString,100);
    % fprintf('Time %s corresponds to frame %d\n', timeString, frameNumber);
    
    % Extract minutes, seconds, and milliseconds
    minutes = str2double(timeString(1:2));
    seconds = str2double(timeString(3:4));
    milliseconds = str2double(timeString(5:6)) * 10;
    
    % Convert to total seconds
    totalSeconds = minutes * 60 + seconds + milliseconds / 1000;
    
    % Convert to frames (FR = 100 fps) and round
    frames = round(totalSeconds * FR);
end

