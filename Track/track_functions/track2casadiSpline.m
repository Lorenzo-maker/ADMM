function track = track2casadiSpline(track)
%
%
%

% create dummy symbolic variable to build functions
alpha_sym = casadi.MX.sym('alpha');

% low-level casadi bspline function call
x_interp = casadi.Function.bspline('x', {track.u}, track.Px, {track.p});
y_interp = casadi.Function.bspline('x', {track.u}, track.Py, {track.p});
z_interp = casadi.Function.bspline('x', {track.u}, track.Pz, {track.p});
wl_interp = casadi.Function.bspline('x', {track.u}, track.Pwl, {track.p});
wr_interp = casadi.Function.bspline('x', {track.u}, track.Pwr, {track.p});
tw_interp = casadi.Function.bspline('x', {track.u}, track.Ptw, {track.p});

% position expression
P = vertcat(x_interp(alpha_sym),y_interp(alpha_sym),z_interp(alpha_sym));
p_interp = vertcat(P , wl_interp(alpha_sym), wr_interp(alpha_sym), tw_interp(alpha_sym));

% tangent vector
T = P.jacobian(alpha_sym);
normts = norm(T);
normts_sq = sum(T.^2);
t = T./normts;

% normal vector
dT = T.jacobian(alpha_sym);
v = vertcat(-t(2), t(1), 0);
vh = v./norm(v);
n = (dT-t*dot(t, dT))./normts_sq;

%curvature
kp = dot(n, vh);

% binormal vector
wh = cross(t, vh);

% ribbon frame
Rgh = [t, vh, wh];
cmu = cos(p_interp(6));
smu = sin(p_interp(6));

% banking transform
Rhs = [1,0,0;0,cmu,-smu;0,smu,cmu];

% tangent frame
Rgs = Rgh*Rhs;

% normal
ns = Rgs(:, 2);

% binormal
ms = Rgs(:, 3);

% ground track transform
g_gs = [Rgs, P; 0,0,0,1];

% normal in ns dir
ns_normal = (dT-t*dot(t, dT))./normts_sq;
kns = dot(ns_normal, Rgs(:, 2));

% angular velocity
b = dot(ns_normal, vh);
c = dot(ns_normal, wh);
omega_x = b*t(3)./(sum(t.^2));
omega_y = -c;
omega_z = b;
omega = vertcat(omega_x, omega_y, omega_z);
omega_gh_h = hat(omega);
mu_dot = jacobian(p_interp(6), alpha_sym)./normts;
omega_hs_h = [0, 0, 0;0, 0, -mu_dot;0, mu_dot, 0];
omega_gs_s = Rhs'*(omega_gh_h + omega_hs_h)*Rhs;
omega_gs_sv = vecForm(omega_gs_s);

% linear velocity
ts_s = Rgs'*t;

% track-induced twist
Vgs_s = [ts_s; omega_gs_sv];



% build functions from expressions
p_fun = casadi.Function('posi', {alpha_sym}, {p_interp}); % position
normts_fun = casadi.Function('normts', {alpha_sym}, {normts}); % norm of tangent vector
t_fun = casadi.Function('t_fun', {alpha_sym}, {t});
Vgs_s_fun = casadi.Function('Vgs_s', {alpha_sym}, {Vgs_s});
mu_dot_fun = casadi.Function('mu_dot', {alpha_sym}, {mu_dot});
kns_fun = casadi.Function('kns', {alpha_sym}, {kns});
ns_fun = casadi.Function('ns', {alpha_sym}, {ns});
ms_fun = casadi.Function('ms', {alpha_sym}, {ms});
Rgs_fun = casadi.Function('Rgs', {alpha_sym}, {Rgs});
g_gs_fun = casadi.Function('g_gs', {alpha_sym}, {g_gs});

% replace all track functions with casadi-bspline-based ones
track.fun_pos = p_fun;
track.fun_g_gs = g_gs_fun;
track.fun_kns = kns_fun;
track.fun_normts = normts_fun;
track.fun_ts = t_fun;
track.fun_Vgs_s = Vgs_s_fun;
track.fun_mudot = mu_dot_fun;
track.fun_ns = ns_fun;
track.fun_ms = ms_fun;
track.fun_Rgs = Rgs_fun;

end