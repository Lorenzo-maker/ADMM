% Import library and global variable
import casadi.*
global car

% Build track
% track = track_fun(1, '3D', 1, 600, 100); % 1-> Nurburgring (600, 100ctrl)
tol = 10^-5; %tollerance to avoid Nan


%% init vehicle model

% Initialize car model
q0 = [0;0;qm;0;0];                                % [e_p0; e_psi0; d0; theta0; phi0]
car = track_vechicle_model(pista, q0);
car.cz1 = cz1;
car.cz2 = cz2;
car.cx = cx;
car.Sf = S;
car.rho = rho;
car.a1 = a1;
car.a2 = a2;
car.q1 = q1;
car.q2 = q2;
car.qm = (car.a2*car.q1 + car.a1*car.q2)/(car.a1+car.a2);
car.h = hg;
car.GM = car.h - car.qm;

% Build inertial and structural car matrix
k1_front = 36000;
k2_rear = 24000;
car.k_phi1 = k_phi1;
car.k_phi2 = k_phi2;
car.k_phi  = car.k_phi1 + car.k_phi2;
car.M = [eye(3).*240, zeros(3,3); zeros(3,3), diag([40, 100, 110])]; % inertia matrix
car.M = adjointStar(RpTohomogeneous(eye(3), [0; 0; GM]))*car.M*adjointInv(RpTohomogeneous(eye(3), [0; 0; GM]));
car.K = [(36000+24000)*2;(36000*2)*a1^2 + (24000*2)*a2^2;car.k_phi];
car.C = [(3280*2 + 2195*2);(3280*2)*a1^2 + (2195*2)*a2^2;0.5*(3280+2195)*t^2];
car.Maxle = zeros(6,6);

car.computeCasadi_state_space_function_alpha_num(); %'compile', true

%% numerical values of car model

g_gs_num = full(pista.fun_g_gs(alfarange));
g_gs_num = reshape(g_gs_num, 4, 4, []);

Jac_num = full(car.T_jac(alfarange));
Jac_der_num = chop(full(car.T_jac_der(alfarange)),1);

%% Definition of states in CasAdi syntax

% symbolic numerical parameters
g_gs_sym = SX.sym('ggs', 4, 4);
Jac_sym = SX.sym('J', 6, 1);
Jac_der_sym = SX.sym('Jdot', 6, 1);
nx = 11; 
X = SX.sym('X', nx);

% curve parameter
% alfa      = X(1);               
% alfa_lb   = 0.;
% alfa_ub   = 1; 

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

%% Definition of inputs in CasAdi syntax" 
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

%% Definition of algebraic parameters in CasAdi syntax" 
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

%% Car's dynamic (system of ODEs)"
alpha_sym = SX.sym('alpha', 1, 1);
qdot  = car.xdot_sym(g_gs_sym, Jac_sym, Jac_der_sym, X, U(1:end-1));

node  = nx;
X_dot = qdot./(qdot(1)+1e-10);
opt = struct;
opt.cse = true;
F = Function('F', {g_gs_sym, Jac_sym, Jac_der_sym, X, U}, {X_dot}); 

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

%% Implementing tyre model Fx

Fzij  = SX.sym('Fzij', 1);
Fz0   = 663.95;             % Nominal wheel load
PDX1  = -2.55009290E+00;  % Lateral friction Muy
PDX2  =  2.45856680E-01;  % Variation of friction Muy with load
% PCX1  =  1.21286210E+00;  % Shape factor Cfy for lateral forces
% PKX1  =  6.61326370E+01;  % Maximum value of stiffness Kfy/Fznom
% PKX2  =  6.26694500E-06;
% PKX3  = -1.80505170E-01;
% PEX1  =  6.38681620E-01;   % Lateral curvature Efy at Fznom
% PEX2  = -1.09337190E-01;   % Variation of curvature Efy with load
% PEX3  = -4.66529940E-01;   % Zero order camber dependency of curvature Efy
% PEX4  = -5.65192410E-02;
% 
% Fzij  = SX.sym('Fzij', 1);
% kx    = SX.sym('kx', 1);
dfz   = (Fzij - Fz0)/Fz0; % Deviation from the nominal tire load
Dx    = (PDX1 + PDX2*dfz)*Fzij;
% Cx    = PCX1;
% Bx    = Fz*(PKX1 + PKX2*dfz)*exp(PKX3*dfz)/(Cx*Dx + tol);
% Ex    = (PEX1 + PEX2*dfz + PEX3*dfz^2)*(1 - PEX4*sign(kx));
% Fxij  = -Dx*sin(Cx*atan(Bx*kx - Ex*(Bx*kx - atan(Bx*kx))));
% 
% Fxtyre = Function('Fxtyre', {kx, Fzij}, {Fxij}, {'kx', 'Fzij'}, {'Fxtyre'}); 
Fxmaxij = Function('Fxmaxij', {Fzij}, {Dx}, {'Fzij'}, {'Dx'}); 

%% Implementing tyre model Fy

Fz0   = 663.95;             % Nominal wheel load
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

%% Load Transfer


% Dz_long = My_actual/l/2;
% Dz_lat  = Mx_actual/t;

%% Slip angle functions 

% Vaxle = car.axle_twist(g_gs_sym, Jac_sym, Jac_der_sym, X);
Vaxle = SX.sym('Vax', 6, 1);
alfa1 = U(7) - (Vaxle(2)+Vaxle(end)*a1)/(Vaxle(1)+tol);
alfa2 = -(Vaxle(2)-Vaxle(end)*a2)/(Vaxle(1)+tol);
Caf = 10e3;
Car = 5e3;
alfa1_fun = Function('alfa1', {Vaxle, U}, {alfa1}); 
alfa2_fun = Function('alfa2', {Vaxle, U}, {alfa2}); 

% Fyn = Fytyre(alfa1_fun(X,U), Fz11tyre(Fz_actual, Dz_long, Dz_lat)) + Fytyre(alfa1_fun(X,U), Fz12tyre(Fz_actual, Dz_long, Dz_lat)) + Fytyre(alfa2_fun(X,U), Fz21tyre(Fz_actual, Dz_long, Dz_lat)) + Fytyre(alfa2_fun(X,U), Fz22tyre(Fz_actual, Dz_long, Dz_lat));
% Mzn = (Fytyre(alfa1_fun(X,U), Fz11tyre(Fz_actual, Dz_long, Dz_lat)) + Fytyre(alfa1_fun(X,U), Fz12tyre(Fz_actual, Dz_long, Dz_lat)))*a1 - (Fytyre(alfa2_fun(X,U), Fz21tyre(Fz_actual, Dz_long, Dz_lat)) + Fytyre(alfa2_fun(X,U), Fz22tyre(Fz_actual, Dz_long, Dz_lat)))*a2;

% Eqns = [Fyn; Mzn];
% Eq = Function('Eq', {X, U}, {Eqns}, {'X', 'U'}, {'Eqns'}); 

%% Lateral load trasnfer function

% Y1 = SX.sym('Y1');
% Y2 = SX.sym('Y2');
% Mx_s = SX.sym('Mx_s');
% a1b = car.a1 - (Y1*a1-Y2*a2)/(Y1+Y2);
% a2b = car.a2 + (Y1*a1-Y2*a2)/(Y1+Y2);
% qb = (a2b*car.q1 + a1b*car.q2)/(a1b+a2b);
% Dz_2 = ((Y1+Y2)*(car.h-qb)*(car.k_phi2/car.k_phi) + Y2*car.q2)/t;
% Dz_1 = ((Y1+Y2)*(car.h-qb)*(car.k_phi1/car.k_phi) + Y1*car.q1)/t;
% DZ_2 = Function('DZ_2', {Y1, Y2}, {Dz_2}, {'Y1', 'Y2'}, {'Dz_2'});
% DZ_1 = Function('DZ_1', {Y1, Y2}, {Dz_1}, {'Y1', 'Y2'}, {'Dz_1'});


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


% Funzione segno
arg = SX.sym('arg');
pos_val = ((tanh(1000*arg))+1)/2;
neg_val = (tanh(1000*arg)-1)/2;
sgn_val = ((tanh(10*arg)));
plus_fun = Function('plus_fun', {arg}, {pos_val});
minus_fun = Function('minus_fun', {arg}, {neg_val});
sgn_fun = Function('sgn_fun', {arg}, {sgn_val});

% Evaluate numerical quantities
tau_colloc = tau_root(2:end);
g_gs_colloc = cell(length(alfarange)-1, length(tau_colloc));
Jac_colloc = cell(length(alfarange)-1, length(tau_colloc));
Jac_der_colloc = cell(length(alfarange)-1, length(tau_colloc));
alfa_colloc = zeros(length(alfarange)-1, d);
for k = 1:length(alfarange)-1
    for j = 1:d
        alfa_colloc(k,j) = min(alfarange(k) + dalfa*tau_colloc(j),1);
        g_gs_colloc{k, j} = full(pista.fun_g_gs(min(alfarange(k) + dalfa*tau_colloc(j),1)));
        Jac_colloc{k, j} = full(car.T_jac(min(alfarange(k) + dalfa*tau_colloc(j),1)));
        Jac_der_colloc{k, j} = chop(full(car.T_jac_der(min(alfarange(k) + dalfa*tau_colloc(j),1))),1);
    end
end
