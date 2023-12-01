function car = car_parameters_ocp_fun(alfarange, pista, colloc_type, d)

import casadi.*
%% Vehicle's data (some of theme are vheicle properties, others are needed to build certain properties)

rho         = 1.2073;                    % air density                (property)
S           = 1.4;                       % front area                 (property)
cx          = 0.84;                      % drag coefficient           (property)
cl          = 1.34;                      % lift coefficient           (property)
ab          = 0.4;                       % aerobalance
G           = 9.81;                      % gravity
m           = 240;                       % mass of car in kg
mu_x        = 1.0 ;                      % dry road example
mu_y        = 1.0;                       % dry road example           (property.data)
tol         = 10^-5;                     % tolerance for avoiding NaN (property.data)
Vi          = 15;                        % initial speed              (property.data)
Vmax        = 150/3.6;                   % max speed                  (property.data)
Pmin        = -1800e3;                   % max braking power          (property.data)
Pmax        = 47e3;                      % max power                  (property.data) 
a1          = 0.765;                                                  %(property)              
a2          = 0.815;                                                  %(property)
cz1         = ab*cl;                                                  %(property)                      
cz2         = (1-ab)*cl;                                              %(property)
t1          = 1.210;
t2          = 1.210;
t           = (t1+t1)/2;                                              %(property.data)
l           = a1 + a2;
q1          = 0.335;                                                  %(property)
q2          = 0.335;                                                  %(property)  
qm          = (a2*q1 + a1*q2)/(a1+a2);                                %(property)
hg          = 0.335+0.1;                                              %(property)
GM          = hg - qm;                                                %(property)
k_phi1      = 0.5*36000*t^2;                                          %(property)  
k_phi2      = 0.5*24000*t^2;                                          %(property)
k1          = 36000;
k2          = 24000;
c1          = 3280;
c2          = 2195;
kb          = 0.6;                                                    %(property.data)

%% init vehicle model

% Initialize car model
q0 = [0;0;qm;0;0];                                % [e_p0; e_psi0; d0; theta0; phi0]
car = track_vechicle_model(pista, q0);
% aero data
car.cz1 = cz1;
car.cz2 = cz2;
car.cx = cx;
car.Sf = S;
car.rho = rho;
% geometric data
car.a1 = a1;
car.a2 = a2;
car.q1 = q1;
car.q2 = q2;
car.qm = (car.a2*car.q1 + car.a1*car.q2)/(car.a1+car.a2);
car.h = hg;
car.GM = car.h - car.qm;
car.data.t = t;
% Powertrain & brake data
car.data.Pmin = Pmin;
car.data.Pmax = Pmax;
car.data.kb = kb;
% Collocation data
car.data.colloc_type = colloc_type;
car.data.d = d;
% Others
car.data.Vi = Vi;
car.data.Vmax = Vmax;
car.data.tol = tol;
% Build inertial car matrix
car.M = [eye(3).*m, zeros(3,3); zeros(3,3), diag([40, 100, 110])]; % inertia matrix
car.M = adjointStar(RpTohomogeneous(eye(3), [0; 0; GM]))*car.M*adjointInv(RpTohomogeneous(eye(3), [0; 0; GM]));
car.Maxle = zeros(6,6);
% Build Structural car matrix
car.k_phi1 = k_phi1;
car.k_phi2 = k_phi2;
car.k_phi  = car.k_phi1 + car.k_phi2;
car.K = [(k1+k2)*2;(k1*2)*a1^2 + (k2*2)*a2^2;car.k_phi];
car.C = [(c1*2 + c2*2);(c1*2)*a1^2 + (c2*2)*a2^2;0.5*(c1+c2)*t^2];

car.computeCasadi_state_space_function_alpha_num(); %'compile', true

%% numerical values of car model

g_gs_num = full(pista.fun_g_gs(alfarange));
g_gs_num = reshape(g_gs_num, 4, 4, []);
car.data.g_gs_num = g_gs_num;

Jac_num = full(car.T_jac(alfarange));
Jac_der_num = chop(full(car.T_jac_der(alfarange)),1);
car.data.Jac_num = Jac_num;
car.data.Jac_der_num = Jac_der_num;


%% Definition of states in CasAdi syntax & scaling factors

% Scaling for states
ep_scale    = 10;
ef_scale    = pi/6;  
d_scale     = 0.1;  
theta_scale = 10*3.14/180;
phi_scale   = 10*3.14/180;
alfa_dot_scale  = 0.0104;
ep_dot_scale    = ep_scale/1;   
ef_dot_scale    = ef_scale/1;  
d_dot_scale     = d_scale/0.1; 
theta_dot_scale = theta_scale/1; 
phi_dot_scale   = phi_scale/1; 

% Symbolic states
nx = 11; 
X = SX.sym('X', nx);

% lateral displacement
ep      = X(1);              
ep_lb   = -1; 
ep_ub   =  1; 

% angular displacement
ef      = X(2);            
ef_lb   = -0.5; 
ef_ub   =  0.5; 

% vertical displacement
dd       = X(3);              
d_lb    = -1; 
d_ub    =  1; 

% pitch
theta       = X(4);             
theta_lb    = -1; 
theta_ub    = 1; 

% roll
phi       = X(5);              
phi_lb    = -1; 
phi_ub    =  1; 

% curve parameter
alfa_dot      = X(6);               
alfa_dot_lb   = 0.;
alfa_dot_ub   = 1; 

% lateral displacement
ep_dot      = X(7);              
ep_dot_lb   = -1; 
ep_dot_ub   =  1; 

% angular displacement
ef_dot      = X(8);            
ef_dot_lb   = -1; 
ef_dot_ub   =  1; 

% vertical displacement
d_dot       = X(9);              
d_dot_lb    = -1; 
d_dot_ub    =  1; 

% pitch
theta_dot       = X(10);             
theta_dot_lb    = -1; 
theta_dot_ub    =  1; 

% roll
phi_dot       = X(11);              
phi_dot_lb    = -1; 
phi_dot_ub    =  1; 

% build vector of boundary and scaling factors
X_lb    = [ep_lb; ef_lb; d_lb; theta_lb; phi_lb; alfa_dot_lb; ep_dot_lb; ef_dot_lb; d_dot_lb; theta_dot_lb; phi_dot_lb];
X_ub    = [ep_ub; ef_ub; d_ub; theta_ub; phi_ub; alfa_dot_ub; ep_dot_ub; ef_dot_ub; d_dot_ub; theta_dot_ub; phi_dot_ub];
X_scale = [ep_scale; ef_scale; d_scale; theta_scale; phi_scale; alfa_dot_scale; ep_dot_scale; ef_dot_scale; d_dot_scale; theta_dot_scale; phi_dot_scale];

% Saving in property.data
car.data.nx = nx;
car.data.X_lb = X_lb;
car.data.X_ub = X_ub;
car.data.X_scale = X_scale;
car.data.ep_scale = ep_scale; 

%% Definition of inputs in CasAdi syntax & scaling factors" 

% Scaling for inputs
Fz_scale = 2*(m*G + 0.5*rho*S*cl*Vmax^2);
Fx_scale = Fz_scale*mu_x;
Fy_scale = Fz_scale*mu_y;
Mx_scale = (Fz_scale*((t1+t2)/2)); 
My_scale = (Fz_scale*((l)/2)); 
Mz_scale = (Fy_scale*(l/2));
delta_scale = 15*3.14/180;

%Symbolic states
nu = 7;  
U = SX.sym('U', nu);

% vertical force road to axle
Fz      = U(3);            
Fz_lb   = 0; 
Fz_ub   = 1; 

% longitudinal force road to axle
Fx      = U(1);         
Fx_lb   = -1; 
Fx_ub   =  1; 

% lateral force road to axle
Fy      = U(2);            
Fy_lb   = -1; 
Fy_ub   =  1; 

% Torque x road to axle
Mx      = U(4);            
Mx_lb   = -1; 
Mx_ub   =  1; 

% Torque y road to axle
My      = U(5);            
My_lb   = -1;
My_ub   =  1; 

% Torque z road to axle
Mz      = U(6);            
Mz_lb   = -1; 
Mz_ub   =  1; 

% Steer angle 
delta      = U(7);            
delta_lb   = -1; 
delta_ub   =  1; 

% build vectors of boundary and scaling factors
U_lb    = [Fx_lb; Fy_lb; 0*Fz_lb; 0*Mx_lb; 0*My_lb; Mz_lb; delta_lb];
U_ub    = [Fx_ub; Fy_ub; 0*Fz_ub; 0*Mx_ub; 0*My_ub; Mz_ub; delta_ub];
U_scale = [Fx_scale; Fy_scale; Fz_scale; Mx_scale; My_scale; Mz_scale; delta_scale];

% Saving in property.data
car.data.nu = nu;
car.data.U_lb = U_lb;
car.data.U_ub = U_ub;
car.data.U_scale = U_scale; 
car.data.Fy_scale = Fy_scale; 
car.data.Mz_scale = Mz_scale; 

%% Definition of algebraic parameters in CasAdi syntax & scaling factors (the same of inputs)" 

%symbolic algebraic
nz = 8;  
Z = SX.sym('U', nz);

% Fz11 
Fz11     = Z(1);            
Fz11_lb  = tol;
Fz11_ub  =  1; 

% Fz12 
Fz12     = Z(2);            
Fz12_lb  = tol;
Fz12_ub  =  1; 

% Fz21 
Fz21     = Z(3);            
Fz21_lb  = tol;
Fz21_ub  =  1; 

% Fz22 
Fz22     = Z(4);            
Fz22_lb  = tol;
Fz22_ub  =  1; 

% Fz22 
Y1tot     = Z(5);            
Y1tot_lb  = -1;
Y1tot_ub  =  1; 

% Fz22 
Y2tot     = Z(6);            
Y2tot_lb  = -1;
Y2tot_ub  =  1;

% Fz22 
Fxa     = Z(7);            
Fxa_lb  = 0;
Fxa_ub  =  1; 

% Fz22 
Fxb     = Z(8);            
Fxb_lb  = -1;
Fxb_ub  =  0;
% build vectors of boundary and scaling factors
Z_lb    = [Fz11_lb; Fz12_lb; Fz21_lb; Fz22_lb; Y1tot_lb; Y2tot_lb; Fxa_lb; Fxb_lb];
Z_ub    = [Fz11_ub; Fz12_ub; Fz21_ub; Fz22_ub; Y1tot_ub; Y2tot_ub; Fxa_ub; Fxb_ub];
Z_scale = [Fz_scale;Fz_scale;Fz_scale;Fz_scale;Fy_scale; Fy_scale; Fx_scale; Fx_scale];

%Saving in property.data
car.data.nz = nz;
car.data.Z_lb = Z_lb;
car.data.Z_ub = Z_ub;
car.data.Z_scale = Z_scale; 
car.data.Fx_scale = Fx_scale; 

%% Car's dynamic (system of ODEs)"

% symbolic numerical parameters
g_gs_sym = SX.sym('ggs', 4, 4);
Jac_sym = SX.sym('J', 6, 1);
Jac_der_sym = SX.sym('Jdot', 6, 1);

% Build casadi Function for dynamics
qdot  = car.xdot_sym(g_gs_sym, Jac_sym, Jac_der_sym, X, U(1:end-1));
X_dot = qdot./(qdot(1)+1e-10);
% opt = struct;
% opt.cse = true;
F = Function('F', {g_gs_sym, Jac_sym, Jac_der_sym, X, U}, {X_dot}); 

%Saving in property.fun
car.fun.F = F;

%% Vertical load for each tyre

Fzij0  = SX.sym('Fzij0', 1);
Myij = SX.sym('Myij',1);
Dz_latij = SX.sym('Dz_latij',1);
Fzaero1 = SX.sym('Fzaero1',1);
Fzaero2 = SX.sym('Fzaero2',1);
Fzaero = Fzaero1 + Fzaero2;
Mzaero = Fzaero2*a2 - Fzaero1*a1;


Fz11_s = (Fzij0-Fzaero)*a2/(a1+a2)/2 + Fzaero1/2 - (Myij - Mzaero)/(2*l) -(Dz_latij);%+k1_front/(k1_front+k2_rear)
Fz11tyre = Function('Fz11tyre', {Fzij0, Fzaero1, Fzaero2, Myij, Dz_latij}, {Fz11_s}, {'Fzij0', 'Fzaero1', 'Fzaero2', 'Myij', 'Dz_latij'}, {'Fz11_s'}); 
Fz12_s = (Fzij0-Fzaero)*a2/(a1+a2)/2 + Fzaero1/2 - (Myij - Mzaero)/(2*l) +(Dz_latij);%k1_front/(k1_front+k2_rear) 
Fz12tyre = Function('Fz12tyre', {Fzij0, Fzaero1, Fzaero2, Myij, Dz_latij}, {Fz12_s}, {'Fzij0', 'Fzaero1', 'Fzaero2', 'Myij', 'Dz_latij'}, {'Fz12_s'}); 
Fz21_s = (Fzij0-Fzaero)*a1/(a1+a2)/2 + Fzaero2/2 + (Myij - Mzaero)/(2*l) -(Dz_latij);%k2_rear/(k1_front+k2_rear)
Fz21tyre = Function('Fz21tyre', {Fzij0, Fzaero1, Fzaero2, Myij, Dz_latij}, {Fz21_s}, {'Fzij0', 'Fzaero1', 'Fzaero2', 'Myij', 'Dz_latij'}, {'Fz21_s'}); 
Fz22_s = (Fzij0-Fzaero)*a1/(a1+a2)/2 + Fzaero2/2 + (Myij - Mzaero)/(2*l) +(Dz_latij);%k2_rear/(k1_front+k2_rear)
Fz22tyre = Function('Fz22tyre', {Fzij0, Fzaero1, Fzaero2, Myij, Dz_latij}, {Fz22_s}, {'Fzij0', 'Fzaero1', 'Fzaero2', 'Myij', 'Dz_latij'}, {'Fz22_s'});  

%Saving in property.fun
car.fun.Fz11tyre = Fz11tyre;
car.fun.Fz12tyre = Fz12tyre;
car.fun.Fz21tyre = Fz21tyre;
car.fun.Fz22tyre = Fz22tyre;


%% Implementing tyre model Fx (Only grip for maximum Fxij)

Fzij  = SX.sym('Fzij', 1);
Fz0   = 663.95;             % Nominal wheel load
PDX1  = -2.55009290E+00;  % Lateral friction Muy
PDX2  =  2.45856680E-01;  % Variation of friction Muy with load
dfz   = (Fzij - Fz0)/Fz0; % Deviation from the nominal tire load
Dx    = (PDX1 + PDX2*dfz)*Fzij;
Fxmaxij = Function('Fxmaxij', {Fzij}, {Dx}, {'Fzij'}, {'Dx'}); 

%Saving in property.fun
car.fun.Fxmaxij = Fxmaxij;

%% Implementing tyre model Fy

Fz0   = 663.95;           % Nominal wheel load
PDY1  = -2.49052010E+00;  % Lateral friction Muy
PDY2  =  1.49411960E-01;  % Variation of friction Muy with load
PCY1  =  1.40000000E+00;  % Shape factor Cfy for lateral forces
PKY1  = -1.04380880E+02;  % Maximum value of stiffness Kfy/Fznom
PKY2  =  3.67230100E+00;
PKY4  =  1;                % Non presente
PEY1  =  2.86078180E-02;   % Lateral curvature Efy at Fznom
PEY2  =  8.58706720E-02;   % Variation of curvature Efy with load
PEY3  = -8.56699590E-01;   % Zero order camber dependency of curvature Efy

alfay = SX.sym('alfay', 1);
dfz   = (Fzij - Fz0)/Fz0; % Deviation from the nominal tire load
Dy    = (PDY1 + PDY2*dfz)*Fzij;
Cy    = PCY1;
By    = Fz0*PKY1*sin(2*atan(Fzij/(Fz0*PKY2)))/(Cy*Dy + tol);
Ey    = (PEY1 + PEY2*dfz)*(1 - PEY3*if_else_smooth(alfay,-alfay,1,-1,'C',20));%sign(alfay)
Ey2   = (PEY1 + PEY2*dfz)*(1 - PEY3*sign(alfay));%
Fyij  = -Dy*sin(Cy*atan(By*alfay - Ey*(By*alfay - atan(By*alfay))));
Fyij2 = -Dy*sin(Cy*atan(By*alfay - Ey2*(By*alfay - atan(By*alfay))));

Fytyre = Function('Fytyre', {alfay, Fzij}, {Fyij}, {'alfay', 'Fzij'}, {'Fyij'}); 
Fytyre2 = Function('Fytyre2', {alfay, Fzij}, {Fyij2}, {'alfay', 'Fzij'}, {'Fyij2'}); 
Fymaxij = Function('Fymaxij', {Fzij}, {Dy}, {'Fzij'}, {'Dy'}); 
Cs = Fz0*PKY1*sin(2*atan(Fzij/(Fz0*PKY2)));
CsFz =  Function('Fymaxij', {Fzij}, {Cs}, {'Fzij'}, {'Cs'});

%Saving in property.fun
car.fun.Fytyre = Fytyre;
car.fun.Fytyre2 = Fytyre2;
car.fun.Fymaxij = Fymaxij;
car.fun.CsFz = CsFz;


%% Slip angle functions 

Vaxle = SX.sym('Vax', 6, 1);
alfa1 = U(7) - (Vaxle(2)+Vaxle(end)*a1)/(Vaxle(1)+tol);
alfa2 = -(Vaxle(2)-Vaxle(end)*a2)/(Vaxle(1)+tol);
Caf = 10e3;
Car = 5e3;
alfa1_fun = Function('alfa1', {Vaxle, U}, {alfa1}); 
alfa2_fun = Function('alfa2', {Vaxle, U}, {alfa2}); 

%Saving in property.fun
car.fun.alfa1_fun = alfa1_fun;
car.fun.alfa2_fun = alfa2_fun;

%% Lateral load trasnfer function

Y1 = SX.sym('Y1');
Y2 = SX.sym('Y2');
Mx_s = SX.sym('Mx_s');
Delta_rs = -Mx_s - (Y1*car.q1 + Y2*car.q2);
a1b = car.a1 - (Y1*a1-Y2*a2)/(Y1+Y2);
a2b = car.a2 + (Y1*a1-Y2*a2)/(Y1+Y2);
qb = (a2b*car.q1 + a1b*car.q2)/(a1b+a2b);
Dz_2 = (Delta_rs*(car.k_phi2/car.k_phi) + Y2*car.q2)/t;
Dz_1 = (Delta_rs*(car.k_phi1/car.k_phi) + Y1*car.q1)/t;
DZ_2 = Function('DZ_2', {Mx_s, Y1, Y2}, {Dz_2}, {'Mx_s', 'Y1', 'Y2'}, {'Dz_2'});
DZ_1 = Function('DZ_1', {Mx_s, Y1, Y2}, {Dz_1}, {'Mx_s', 'Y1', 'Y2'}, {'Dz_1'});

%Saving in property.fun
car.fun.DZ_1 = DZ_1;
car.fun.DZ_2 = DZ_2;

%% Evaluate numerical quantities & save in property.data (maybe they can be evaluated in sub_opti_map ?)
tau_root = [0 collocation_points(d, colloc_type)];
tau_colloc = tau_root(2:end);

pre = diff(alfarange);
index = find(pre<0);
if ~isempty(index)
    pre(index) = 1-alfarange(index);
end
dalfa_vec = pre;
car.data.dalfa_vec = dalfa_vec;   
g_gs_colloc = cell(length(alfarange)-1, length(tau_colloc));
Jac_colloc = cell(length(alfarange)-1, length(tau_colloc));
Jac_der_colloc = cell(length(alfarange)-1, length(tau_colloc));
alfa_colloc = zeros(length(alfarange)-1, d);
for k = 1:length(alfarange)-1
    for j = 1:d
        alfa_colloc(k,j) = min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1);
        g_gs_colloc{k, j} = full(pista.fun_g_gs(min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1)));
        Jac_colloc{k, j} = full(car.T_jac(min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1)));
        Jac_der_colloc{k, j} = chop(full(car.T_jac_der(min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1))),1);
    end
end
car.data.g_gs_colloc = g_gs_colloc;
car.data.Jac_colloc = Jac_colloc;
car.data.Jac_der_colloc = Jac_der_colloc;

end

