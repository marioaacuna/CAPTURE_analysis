function vel = calculate_velocity(position, fps)
    % Calculate velocity using central difference
    vel = diff(position) * fps;
    % Add duplicate of last velocity to match original length
    vel = [vel; vel(end,:)];
end