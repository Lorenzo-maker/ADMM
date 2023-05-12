function alpha_new = trackCurvatureAdaptiveMesh(track, delta_len_ref, n_skip_straight, n_skip_curve, curve_buffer, options)
    arguments
        track; delta_len_ref; n_skip_straight; n_skip_curve; curve_buffer;
        options.alpha_in = 0; options.alpha_end = 1; options.curv_treshold = [];
    end
% - OPTIONS: 
%       'alpha_in' = starting alpha
%       'alpha_end' = ending alpha
%       'curv_treshold' = Custom value for treshold curvature
% - track = track model obtained by nurbs_ribbon_casadi
% - delta_len_ref = distance (meters) between two consecutive point on the grid
% - n_skip_straight = number of mesh point skipped in straight (if delta_len_ref = 5, and n_skip_straight = 2, you have one point every 10 meters in the straights) 
% - n_skip_curve = number of mesh point skipped in curve
% - curve_buffer = define the distance before and after a curve fro which we
%   want the same mesh of the correspoding curve (i.e. if delta_len_ref = 5
%   and curve_buffer = 3 in the 15 meters before and after the curve, we have
%   the same mesh grid of the turn


alpha_vec = linspace(options.alpha_in, options.alpha_end, 1500);
L = full(track.length_eval(alpha_vec));
alpha_vec = linspace(options.alpha_in, options.alpha_end, round(L/delta_len_ref));


k = abs(full(track.fun_kp(alpha_vec)));
if isempty(options.curv_treshold)
    curv_treshold = mean(k)*1.2;
else
    curv_treshold = options.curv_treshold;
end
id = k > curv_treshold;



alpha_new = [];

straightOrCurve = k(1)> curv_treshold;  % if true, we start on a curve


id_start = 1;
while 1

    if ~straightOrCurve 
        index = find(id, 1);
        id_end = max(index - curve_buffer, 1);
        if id_end > length(alpha_vec)
            id_end = length(alpha_vec);
        end

        if isempty(index)
            id_end = length(alpha_vec);
            alpha_new = [alpha_new, alpha_vec(id_start:n_skip_straight+1:id_end)];
            if alpha_new(end) ~= options.alpha_end
                alpha_new = [alpha_new, options.alpha_end];
            end
            return
        end
        alpha_new = [alpha_new, alpha_vec(id_start:n_skip_straight+1:id_end)];
    else 
        index = find(~id, 1);
        id_end = index + curve_buffer;
        if id_end > length(alpha_vec)
            id_end = length(alpha_vec);
        end
        if isempty(index)
            id_end = length(alpha_vec);
            alpha_new = [alpha_new, alpha_vec(id_start:n_skip_curve+1:id_end)];
            if alpha_new(end) ~= options.alpha_end
                alpha_new = [alpha_new, options.alpha_end];
            end
            return
        end
        alpha_new = [alpha_new, alpha_vec(id_start:n_skip_curve+1:id_end)];
    end
%     id_start = id_end + 1;
    id = id(id_end+1:end);
    alpha_vec = alpha_vec(id_end+1:end);
    straightOrCurve = ~straightOrCurve;

end

end