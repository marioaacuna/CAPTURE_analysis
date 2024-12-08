function logger(message, log_type)
    if nargin < 2
        log_type = 'INFO';
    end
    
    % Get the current time and user information
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    user = getenv('USER');
    if isempty(user)
        user = getenv('USERNAME');
    end
    computer_name = getenv('COMPUTERNAME');
    if isempty(computer_name)
        computer_name = getenv('HOSTNAME');
    end
    
    % Create the log message
    log_message = sprintf('[%s] [%s] [%s@%s] %s\n', timestamp, log_type, user, computer_name, message);
    
    % Append the log message to the log file
    log_file = fullfile(fileparts(mfilename('fullpath')), 'analysis_log.txt');
    fid = fopen(log_file, 'a');
    if fid == -1
        error('Cannot open log file: %s', log_file);
    end
    fprintf(fid, '%s', log_message);
    fclose(fid);
    
    % Display the log message in the command window
    fprintf('%s', log_message);
end