function psi_track = atan_track_diff(ts_track, direction, offset)
    arguments
        ts_track; direction = 'clockwise'; offset = 0
    end
    ts_track = ts_track(1:2,:);
    norm_ts = vecnorm(ts_track);
    ts_track = ts_track./norm_ts;
    offset_track = atan2(ts_track(2,1),ts_track(1,1));
    ns_track = [-ts_track(2,:);ts_track(1,:)];
    diff_ts = diff(ts_track,1,2);
    sign_dpsi = sign(dot(diff_ts, ns_track(:,2:end)));
    delta_psi = sign_dpsi.*acos(dot(ts_track(:,1:end-1), ts_track(:,2:end)));

    if offset == 0
        psi_track = [offset_track + offset, offset_track + offset + cumsum(delta_psi)];
    else
        psi_track = [offset, offset + cumsum(delta_psi)];
    end


end