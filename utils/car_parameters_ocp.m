
%% Dati
rho         = 1.2073;                    % air density
S           = 1.4;                       % front area
cx          = 0.84;%1.19;                      % drag coefficient
cl          = 1.34;%1.86;                      % lift coefficient
ab          = 0.4;
G           = 9.81;                      % gravity
m           = 240;%647.7;                     % mass of car in kg
Fz_start    = m*G;                       % weigth force (Newton)
mu_x        = 1.0 ;                      % dry road example
mu_y        = 1.0;                       % dry road example
tol         = 10^-3;                     % tolerance for avoiding NaN
Vm          = 22.23;                     % mean total speed (m/s) ~ 80 km/h
hm          = 0.25;                      % mean discretization 
Vi          = 15;                        % initial speed
Vmax        = 150/3.6;                   % max speed
Pmin        = -1800e3;                   % max braking power
Pmax        = 47e3;                      % max power
a1          = 0.765;
a2          = 0.815;
cz1         = ab*cl;
cz2         = (1-ab)*cl;
t1          = 1.210;
t2          = 1.210;
t           = (t1+t1)/2;
l           = a1 + a2;
q1          = 0.335;
q2          = 0.335;
qm          = (a2*q1 + a1*q2)/(a1+a2);
hg          = 0.335+0.1;
GM          = hg - qm;
k_phi1      = 0.5*36000*t^2;
k_phi2      = 0.5*24000*t^2;
k_phi       = k_phi1 + k_phi2;
kb          = 0.6;

%% Scaling factors
%ds = full(pista.length_eval(linspace(0,dalfa,1000)));

% q
alfa_scale  = 1; %1;
ep_scale    = 10;
ef_scale    = pi/6; %1; 
d_scale     = 0.1; %0.5; 
theta_scale = 10*3.14/180;
phi_scale   = 10*3.14/180;

% qdot
alfa_dot_scale  = 0.0104;%dalfa/(4/(1.5*Vmax)); %1;
ep_dot_scale    = ep_scale/1; %10/ep_scale/0.01;  
ef_dot_scale    = ef_scale/1;  %pi/2/ef_scale/0.01;
d_dot_scale     = d_scale/0.1; %0.1/0.01;
theta_dot_scale = theta_scale/1; %10*3.14/180/0.01;
phi_dot_scale   = phi_scale/1; %10*3.14/180/0.01;

X_scale = [ep_scale; ef_scale; d_scale; theta_scale; phi_scale; alfa_dot_scale; ep_dot_scale; ef_dot_scale; d_dot_scale; theta_dot_scale; phi_dot_scale];
nx = length(X_scale);

% if alpha_numeric 
%     nx = length(X_scale)-1; %peculiarity of this model that has alfa numerical
%     nx_full = length(X_scale);
% else
%     nx = length(X_scale);
%     nx_full = length(X_scale);
% end

% Input
Fz_scale = 2*(m*G + 0.5*rho*S*cl*Vmax^2);
Fx_scale = Fz_scale*mu_x;
Fy_scale = Fz_scale*mu_y;
Mx_scale = (Fz_scale*((t1+t2)/2)); 
My_scale = (Fz_scale*((l)/2)); 
Mz_scale = (Fy_scale*(l/2));
delta_scale = 15*3.14/180;

% Others
r_scale  = 4;
u_scale  = Vmax;
v_scale  = 6;
P_scale  = abs(Pmin);

U_scale = [Fx_scale; Fy_scale; Fz_scale; Mx_scale; My_scale; Mz_scale; delta_scale];
Z_scale = [Fz_scale;Fz_scale;Fz_scale;Fz_scale;Fy_scale; Fy_scale; Fx_scale; Fx_scale];
nu = length(U_scale);
nz = length(Z_scale);

