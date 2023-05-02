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
