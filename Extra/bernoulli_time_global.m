import casadi.*
clear all
clc
%% Dati problema 

G = 9.81;
tol = 1e-3;

%% Griglia

N = 100;
d = 2;
p = [];%MX.sym('p');
%% Stati simbolici (nx, X, X_lb, X_ub, Xb) 
nx = 3;
X = MX.sym('X', nx);

% Posizione x
x = X(1);
x_lb = 0;
x_ub = 2;

% Posizione y
y = X(2);
y_lb = -1;
y_ub = 1;

% velocità teta
v = X(3);
v_lb = 0;
v_ub = 10;

% Vettori lower-upper bounds
X_lb = [x_lb; y_lb; v_lb];
X_ub = [x_ub; y_ub; v_ub];

% Struttura per NLP
Xb.lb = X_lb;
Xb.ub = X_ub;

%% Controlli Simbolici (nu, U, U_lb, U_ub, Ub)
nu = 1;
U = MX.sym('U', nu);

teta = U(1);
teta_lb = 0;
teta_ub = pi;

U_lb = teta_lb;
U_ub = teta_ub;

% Struttura per NLP
Ub.lb = U_lb;
Ub.ub = U_ub;

%% Parametri algebrici simbolici (nz, Z, Z_lb, Z_ub, Zb)

nz = 1;
Z = MX.sym('Z', nz);

% dt
h = Z(1);
h_lb = tol;
h_ub = 1;

Z_lb = h_lb;
Z_ub = h_ub;

% Struttura per NLP
Zb.lb = Z_lb;
Zb.ub = Z_ub;

%% Dinamica e Funzione campo (F)

x_dot = v*sin(teta);
y_dot = v*cos(teta);
v_dot = G*cos(teta);

X_dot = [x_dot; y_dot; v_dot];
node  = 3;

% Funzione casadi per dinamica

F = Function('F', {X, U}, {X_dot}, {'X', 'U'}, {'X_dot'}); 

%% Funzione di costo 
Lj = Z^2;
L = Function('L', {Z}, {Lj}, {'Z'}, {'Lj'}); 

%% Condizioni iniziali (C.I. @ t = 0) (X_0)

% states
x_0     = 0;
y_0     = 0;
v_0     = 0;

X_0     = [x_0; y_0; v_0];

% Struttura per NLP
Xb.x0 = X_0;
%% Condizioni finali (C.F. @ t= tf) (X_f)

x_f    = 1;
y_f    = 0.5;

% parametro utile per valutare initial guess
teta_line = atan((x_f - x_0)/(y_f - y_0));


%% Build problem
tic
% Compute w, lbw, ubw and define useful symbolic variables

Pdata.Pg = zeros(1, N);
Pdata.Pgb = [-rand(1, N)*100 ; rand(1, N)*100];
nuc = 0;
nzc = 0;
Wgb.lb = 1e-9;
Wgb.ub = inf;
% pb = nlp_casadi_MX(nx, nu, nz, N, d, Xb, Ub, Zb, 1, 1, Pdata, 1, p);
pb = nlp_casadi_2(nx, nu, 0, nuc, nzc, 1, N, d, Xb, Ub, [], [], [], Wgb, 1, 2, Pdata, 1, 0, 'SX');

% Compute w0 (Define Initial Guess)
k = linspace(1, N, N);
X_0 = k.*[x_f/N; y_f/N; 0] + [0; 0; 1];
U_0 = teta_line*ones(1,N);
Z_0 = (1/N)*ones(1,N);
pb.w0 = [U_0; repmat(X_0, d, 1); repmat(U_0, d, nuc); repmat(U_0, d, nzc); X_0];
pb.w0 = [X_0(:,1); pb.w0(:); 10];

% Compute J, g, lbg, ubg by appending constraints for a generic interval
% Start constraint
% pb.append_g((pb.x_1 - [0;0;0]), zeros(3,1), zeros(3,1), 'start');
% Collocation equations for generic interval
for i = 1:d
   pb.append_g(F(pb.xc(:,i), pb.u )*(pb.wg/N) - pb.xp(:,i), zeros(nx,1), zeros(nx,1)) ;
   pb.Jj = [pb.Jj; L(pb.wg/N)];
end
% Evaluate cost function for generic interval
pb.cost;
% Continuity equation for generic interval
pb.append_g(pb.xc_end - pb.x , zeros(nx,1), zeros(nx,1));
% Equality h for generic interval
% pb.append_g(pb.z_1 - pb.z, 0, 0);
%random constraint
% pb.append_g(pb.x(3) , pb.pgb(1), pb.pgb(2));
% Terminal constraints 
pb.append_g(pb.x(1:2) - [x_f;y_f], zeros(2,1), zeros(2,1), 'end');
% Build map function 
pb.build_map;
% Build g 
pb.build_g;
% Build J
pb.build_J;

% pb = nlp_casadi_SX_2();

% problem construction
prob = struct('f', pb.J, 'x', pb.w, 'g', pb.g);
opts = struct;
opts.ipopt.linear_solver = 'ma57';
opts.ipopt.max_iter = 5000;
opts.ipopt.print_level = 5;
solver = nlpsol('solver', 'ipopt', prob, opts);
toc

% Solve NLP
sol_sym       = solver('x0', pb.w0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg);
siz           = length(sol_sym.x);  
sol.x_opt     = full(sol_sym.x(1:(nx+nu+(nx+nuc+nzc)*d):siz));
sol.y_opt     = full(sol_sym.x(2:(nx+nu+(nx+nuc+nzc)*d):siz));
sol.v_opt     = full(sol_sym.x(3:(nx+nu+(nx+nuc+nzc)*d):siz));
sol.teta_opt  = full(sol_sym.x(4:(nx+nu+(nx+nuc+nzc)*d):siz));
sol.time_opt  = full(sol_sym.x(end)); 
close all
plot(sol.x_opt, -sol.y_opt)
