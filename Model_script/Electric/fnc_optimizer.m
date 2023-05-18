function opt = fnc_optimizer(options,model,sys,rel,driver,resume_sim,resume_out_file)
%
import casadi.*
%
% Objective term
opt.L_kj = ...
    driver.K1 * sys.parameters.dt ...
    + driver.K2 * (rel.der_r / driver.der_r_factor)^2 ...
    + driver.K3 * (rel.der_v / driver.der_v_factor)^2;
%
% Continuous time dynamics
opt.f    = Function('f', {sys.states.X, sys.inputs.U, sys.parameters.Z, sys.gravity.g3D}, {rel.X_dot, opt.L_kj});
opt.Eq   = Function('Eq',{sys.states.X, sys.inputs.U, sys.parameters.Z, sys.gravity.g3D}, {rel.eq}, {'X','U','Z','sys.gravity.g3D'}, {'eq'});
%Function to compute 3D gravity
Rotsg = SX.sym('Rotsg',3,3);
g3D_sym   = [cos(sys.states.psi), sin(sys.states.psi), 0; -sin(sys.states.psi), cos(sys.states.psi), 0; 0, 0, 1]*Rotsg*[0,0,model.global.g]';
fun_g3D = Function('g3D',{sys.states.psi_norm, Rotsg}, {g3D_sym}, {'sys.states.psi_norm','Rotsg'}, {'g3D_sym'});
% Definition of collocation points
col = col_points(model.opt.d_colloc);
%
if resume_sim
    resume_out = load(resume_out_file, 'out');
    resume_out = resume_out.out;
end
%% Build NLP
disp('Building NLP')
% Compute w, lbw, ubw and define useful symbolic variables
opt.Xb.lb = sys.states.X_lb;
opt.Xb.ub = sys.states.X_ub;
opt.Xb.x0 = sys.states.X_0;
opt.Ub.lb = sys.inputs.U_lb;
opt.Ub.ub = sys.inputs.U_ub;
opt.Zb.lb = sys.parameters.Z_lb;
opt.Zb.ub = sys.parameters.Z_ub;
opt.Npg = 3*3 + 2 + 2;
opt.Npgb = 2;
opt.Nthread = 1;
opt.p = 2;

for k = 1:model.opt.dpts
    opt.data_g.Pg(:,k) = [model.track.RVSG{k+1}(:); model.track.pos_grid(1:2,k+1); model.track.ns_grid(1:2,k+1)];
end
opt.data_g.Pgb = [-(model.track.pos_grid(4,2:end)-model.opt.ep_margin);model.track.pos_grid(4,2:end)-model.opt.ep_margin];
% pb = nlp_casadi_SX(sys.states.nx, sys.inputs.nu, sys.parameters.np, model.opt.dpts, model.opt.d_colloc, opt.Xb, opt.Ub, opt.Zb, opt.Npg, opt.Npgb, opt.data_g, opt.Nthread, opt.p);
% pb = nlp_casadi_MX(sys.states.nx, sys.inputs.nu, sys.parameters.np, model.opt.dpts, model.opt.d_colloc, opt.Xb, opt.Ub, opt.Zb, opt.Npg, opt.Npgb, opt.data_g, opt.Nthread, opt.p);
% pb = nlp_casadi_MXSX(sys.states.nx, sys.inputs.nu, sys.parameters.np, model.opt.dpts, model.opt.d_colloc, opt.Xb, opt.Ub, opt.Zb, opt.Npg, opt.Npgb, opt.data_g, opt.Nthread);
pb = nlp_casadi_2(sys.states.nx, sys.inputs.nu, sys.parameters.np, 0, 0, 0, model.opt.dpts, model.opt.d_colloc, opt.Xb, opt.Ub, opt.Zb, [], [], [], opt.Npg, opt.Npgb, opt.data_g, opt.Nthread, opt.p, 'SX');

% Compute w0 (Define Initial Guess)
for k = 1:model.opt.dpts
    initial_guess
    if k == 1
        if isempty(opt.Xb.x0)
            pb.w0 = [pb.w0; sys.states.X_0; sys.inputs.U_0; sys.parameters.Z_0; repmat(sys.states.X_0, model.opt.d_colloc+1, 1)];
        else
            pb.w0 = [pb.w0; opt.Xb.x0; sys.inputs.U_0; sys.parameters.Z_0; repmat(sys.states.X_0, model.opt.d_colloc+1, 1)];
        end
    else
        pb.w0 = [pb.w0; sys.inputs.U_0; sys.parameters.Z_0; repmat(sys.states.X_0, model.opt.d_colloc+1, 1)]; 
    end
end

% Compute J, g, lbg, ubg by appending constraints for a generic interval
% Collocation equations for generic interval
Rgs = reshape(pb.pg(1:9),3,3);
for i = 1:model.opt.d_colloc
    g3DKj = fun_g3D(pb.xc(6,i),Rgs);
    [opt.fr, opt.qr] = opt.f(pb.xc(:,i), pb.u, pb.z, g3DKj);
    pb.append_g(opt.fr*pb.z(end)*sys.parameters.dt_scaling - pb.xp(:,i).*sys.states.X_scale, zeros(pb.nx,1), zeros(pb.nx,1)) ;
    pb.Jj = [pb.Jj; opt.qr * pb.z(end)*sys.parameters.dt_scaling];
end
pb.Jj = pb.Jj*(pb.p(1)) + (1-pb.p(1))*(pb.z(11)*sys.parameters.ep_scaling)^2;
% Evaluate cost function for generic interval
pb.cost;
% Driver
pb.append_g((pb.u(3) - pb.u_1(3))*sys.inputs.delta_v_scaling/(pb.z(end)*sys.parameters.dt_scaling + 1e-5), -pi, pi);
% Continuity equation for generic interval
pb.append_g(pb.xc_end - pb.x, zeros(pb.nx,1), zeros(pb.nx,1));
% Throttle/Braking
pb.append_g(pb.u(1)*pb.u(2), 0, 1e-5);
% beta_r
pb.append_g(pb.z(9) - pb.u(2), -10, 0);
% gamma_r
pb.append_g(pb.z(10) - pb.u(1), -10, 0);
%adherence_f
pb.append_g((pb.z(1) * sys.parameters.Fxf_scaling / ((model.vehicle.tire.mux_0(1) + model.vehicle.tire.k_mux(1) * pb.z(5) * sys.parameters.Fzf_scaling) * pb.z(5) * sys.parameters.Fzf_scaling) ...
        * 0.5 * (model.opt.gripxa_max(1)^-1 + model.opt.gripxb_max(1)^-1 + if_else_smooth(pb.z(1),0,1,-1) * (model.opt.gripxa_max(1)^-1 - model.opt.gripxb_max(1)^-1)))^2 ...
        + (pb.z(3) * sys.parameters.Fyf_scaling / (model.opt.gripy_max(1) * (model.vehicle.tire.muy_0(1) + model.vehicle.tire.k_muy(1) * pb.z(5) * sys.parameters.Fzf_scaling) * pb.z(5) * sys.parameters.Fzf_scaling))^2, 0, 1);
%adherence_r
pb.append_g((pb.z(2) * sys.parameters.Fxr_scaling / ((model.vehicle.tire.mux_0(2) + model.vehicle.tire.k_mux(2) * pb.z(6) * sys.parameters.Fzr_scaling) * pb.z(6) * sys.parameters.Fzr_scaling) ...
        * 0.5 * (model.opt.gripxa_max(2)^-1 + model.opt.gripxb_max(2)^-1 + if_else_smooth(pb.z(2),0,1,-1) * (model.opt.gripxa_max(2)^-1 - model.opt.gripxb_max(2)^-1)))^2 ...
        + (pb.z(4) * sys.parameters.Fyr_scaling / (model.opt.gripy_max(2) * (model.vehicle.tire.muy_0(2) + model.vehicle.tire.k_muy(2) * pb.z(6) * sys.parameters.Fzr_scaling) * pb.z(6) * sys.parameters.Fzr_scaling))^2, 0, 1);
%ep
pb.append_g(pb.z(11)*sys.parameters.ep_scaling, pb.pgb(1), pb.pgb(2));
%track
pos = pb.pg(10:11);
n = pb.pg(12:13);
pb.append_g(pb.x(4:5) - (pos + n*pb.z(11)*sys.parameters.ep_scaling)./sys.states.X_scale(4:5), zeros(2,1), zeros(2,1))
%Eq
g3DK = fun_g3D(pb.x(6),Rgs);
pb.append_g(vertcat(opt.Eq(pb.x, pb.u, pb.z, g3DK)), zeros(rel.neq,1), zeros(rel.neq,1));
%terminal
pb.append_g((model.vehicle.bat.SOC_min + pb.x(7) * sys.states.SOC_scaling) - pb.p(2)*model.opt.SOC_final, 0, 1, 'end');
% % Cyclic
% pb.append_g([pb.x0(1:3);pb.z0(11)]-[pb.x_end(1:3);pb.z_end(11)], [-1;zeros(3,1)], [1;zeros(3,1)], 'cyclic');
% Build map function 
pb.build_map('JN', pb.N);
% Build g 
pb.build_g;
% Build J
pb.build_J('range', 1:pb.N);

opt.d = col.d;

%
%% Create an NLP solver
disp('Starting to build the solver')
% opt.prob = struct('f', pb.J, 'x', pb.w, 'g', pb.g, 'p', pb.p);
% opt.options = options;
% opt.solver = nlpsol('solver', 'ipopt', opt.prob, opt.options );
%
pb.build_solver(options);
opt.solver = pb.solver;
%% Solve the NLP
%
if model.opt.J_Homotopy == 1
    opt.sol = opt.solver('x0', pb.w0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'p', [0;1]);
    opt.sol = opt.solver('x0', opt.sol.x, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'lam_g0', opt.sol.lam_g, 'lam_x0', opt.sol.lam_x, 'p', [0.5,1]);
    opt.sol = opt.solver('x0', opt.sol.x, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'lam_g0', opt.sol.lam_g, 'lam_x0', opt.sol.lam_x, 'p', [1,1]);
elseif model.opt.J_Homotopy == 2
    opt.sol = opt.solver('x0', pb.w0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'p', [1;0]);
    opt.sol = opt.solver('x0', opt.sol.x, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'lam_g0', opt.sol.lam_g, 'lam_x0', opt.sol.lam_x, 'p', [1,0.5]);
    opt.sol = opt.solver('x0', opt.sol.x, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'lam_g0', opt.sol.lam_g, 'lam_x0', opt.sol.lam_x, 'p', [1,0.75]);   
    opt.sol = opt.solver('x0', opt.sol.x, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg',pb.ubg, 'lam_g0', opt.sol.lam_g, 'lam_x0', opt.sol.lam_x, 'p', [1,1]);
else
    if resume_sim
        opt.sol = opt.solver('x0', resume_out.w_opt, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg, 'lam_g0', resume_out.lam_g, 'lam_x0', resume_out.lam_x, 'p', [1;1]);
    else
        opt.sol = opt.solver('x0', pb.w0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg, 'p', [1;1]);
    end
end
opt.w_result = full(opt.sol.x);
opt.lg = full(opt.sol.lam_g);
opt.lx = full(opt.sol.lam_x);
opt.L_result = full(opt.sol.f);
end
%
function o = col_points(d_colloc)
% Collocation points definition
import casadi.*
% Degree of interpolating polynomial
o.d = d_colloc;
% Get collocation points
o.tau_root = [0 collocation_points(o.d, 'legendre')];
% Coefficients of the collocation equation
o.C = zeros(o.d+1,o.d+1);
% Coefficients of the continuity equation
o.D = zeros(o.d+1, 1);
% Coefficients of the quadrature function
o.B = zeros(o.d+1, 1);
% Construct polynomial basis
for j=1:o.d+1
    % Construct Lagrange polynomials to get the polynomial basis at the
    % collocation point
    o.coeff = 1;
    for r=1:o.d+1
        if r ~= j
            o.coeff = conv(o.coeff, [1, -o.tau_root(r)]);
            o.coeff = o.coeff / (o.tau_root(j)-o.tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of
    % the continuity equation
    o.D(j) = polyval(o.coeff, 1.0);
    % Evaluate the time derivative of the polynomial at all collocation
    % points to get the coefficients of the continuity equation
    o.pder = polyder(o.coeff);
    for r=1:o.d+1
        o.C(j,r) = polyval(o.pder, o.tau_root(r));
    end
    % Evaluate the integral of the polynomial to get the coefficients of the
    % quadrature function
    o.pint = polyint(o.coeff);
    o.B(j) = polyval(o.pint, 1.0);
end
%
end
%