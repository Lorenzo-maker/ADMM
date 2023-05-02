function pista = track_fun(index, geometry, start, sector, ctrlpts, NameValueArgs)%degree, approx, eval_points, knots_method)
% track = track_fun(index, geometry, start, sector, ctrlpts, NameValueArgs)
%
% track_fun is a matlab function necessary to preprocess track files and to load them 

% index: enter an index to choose the track
%   - 0 : for Barcelona
%   - 1 : for Nordschleife
%   - 2 : for Calabogie
%    
% geometry: enter the type of geometry track
%   - '2D' : for two-dimensional track
%   - '3D' : for three-dimensional track 
%    
% start: enter the starting point of the track
%
% sector: enter the ending point of the track (enter 'end' to get to the end)
%
% ctrlpts: enter the number of spline control points
%
% NameValueArgs: input with a default value
%   - NameValueArgs.degree : the spline degree (default 3)
%   - NameValueArgs.approx : the method to approximate the spline 
%                               - 'interp' : to interpolate the track with the spline
%                               - 'fit' : to fit the track with the spline (default)
%   - NameValueArgs.eval_points : the method to calculate the point for imposing equalities (Nurbs book pag.364-365)
%                                   - chord (default) 
%                                   - centripetal 
%   - NameValueArgs.knot_metod: the method to calculate the knot points (Nurbs book pag.365/pag.412)
%                                   - equally (available only for interpolating problem)
%                                   - averaging (available only for interpolating problem)
%                                   - deboor (available only for fitting problem) (default)




% preprocessing of track data... we could develop a method for deleting
% some track point in a smart way (Douglas-Peucker)

    arguments
        index;
        geometry;
        start;
        sector;
        ctrlpts;
        NameValueArgs.degree = 3;
        NameValueArgs.approx = 'fit';
        NameValueArgs.eval_points = 'chord';
        NameValueArgs.knot_method = 'deboor';
        NameValueArgs.alfa_limited = false;
        NameValueArgs.alfa_lim = [0, 1];
    end
    % to change these parameters call function(...,'approx','interp',...)
    
    import casadi.*
    if index == 0 
        load SPAn2.mat SPAn2;
        track = SPAn2(1:3:end,:);
%         load Catalunya.mat Catalunya;
%         track = Catalunya(1:2:end,:);
    elseif index == 1
        load Nordschleife.mat drdTable;
        track = [drdTable.x_coord drdTable.y_coord drdTable.z_coord drdTable.width/2 drdTable.width/2 atan(drdTable.lat_inc)];
        track = track(1:4:end,:);
        if geometry == '2D' %#ok<*BDSCA>
            track(:,3) = track(:,3);
            track(:,6) = track(:,6);
        end
    elseif index == 2
        load Calabogie.mat drdTable;
        track = [drdTable.x_coord drdTable.y_coord drdTable.z_coord drdTable.width/2 drdTable.width/2 atan(drdTable.lat_inc)];
        track = track(1:4:end,:);
        if geometry == '2D' 
            track(:,3) = 0*track(:,3);
            track(:,6) = 0*track(:,6);
        end
    elseif index == 3
    load Melbourne.mat Melbourne;
    track = Melbourne(1:3:end,:);
        if geometry == '2D' 
            track(:,3) = 0*track(:,3);
            track(:,6) = 0*track(:,6);
        end
    end
    
    pista = nurbs_ribbon_casadi(track, start, sector, NameValueArgs.degree, ctrlpts, NameValueArgs.approx, NameValueArgs.eval_points, NameValueArgs.knot_method, NameValueArgs.alfa_limited, NameValueArgs.alfa_lim);
    pista.solve();

    %Round control and knots points (useful in order to construct SX NLP and solve numerical issues) 
    p.track.round_factor = 10^14;
    pista.Px  = round(pista.Px.*p.track.round_factor)./p.track.round_factor;
    pista.Py  = round(pista.Py.*p.track.round_factor)./p.track.round_factor;
    pista.Pz  = round(pista.Pz.*p.track.round_factor)./p.track.round_factor;
    pista.Pwl = round(pista.Pwl.*p.track.round_factor)./p.track.round_factor;
    pista.Pwr = round(pista.Pwr.*p.track.round_factor)./p.track.round_factor;
    pista.Ptw = round(pista.Ptw.*p.track.round_factor)./p.track.round_factor;
    pista.u   = round(pista.u.*p.track.round_factor)./p.track.round_factor;


%     pista.eval(SX.sym('us'));
    %pista.rnormal(SX.sym('us'));
    pista.Sframe(SX.sym('us')); 
end