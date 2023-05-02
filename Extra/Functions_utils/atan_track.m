function psi_track = atan_track(ts_track, direction)
    if strcmp(direction, 'clockwise')
        offset = -2*pi;
        offset_diff = -2*pi;
    else 
        offset = 0;
        offset_diff = 2*pi;
    end
    
    for i = 1:length(ts_track(1,:))
        psi_track(i) = mod(atan2(ts_track(2,i),ts_track(1,i)),offset);
        if i > 1
            if abs(psi_track(i) - psi_track(i-1)) > pi
                psi_track(i) = psi_track(i) + offset_diff;
            end
        end
    end


end