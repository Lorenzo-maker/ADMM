clear all
clc
close all
import casadi.*

%% Posture

N = 5;
l = 2;
mb1 = 1;
mb2 = 1;
md = 1;
rd = 0.5;


q0 = [0;0;0;0;0];

q = SX.sym('q', 5, 1);
g0{1} = TrotZ(q0(1));
g0{2} = TrotY(q0(2));
g0{3} = TrotX(q0(3));
g0{4} = Ttz(-l)*TrotX(q0(4));
g0{5} = Tty(l)*TrotY(q0(5));

X = [[0;0;0;0;0;1],...
 [0;0;0;0;1;0],...
 [0;0;0;1;0;0],...
 [0;0;0;1;0;0],...
 [0;0;0;0;1;0]];

[B, Blocal] = FWkin_localPOE(g0, X, q);
B_sym = Function('B_sym', {q}, {B});
B3_sym = Function('B3_sym', {q}, {Blocal{1}*Blocal{2}*Blocal{3}});
B4_sym = Function('B4_sym', {q}, {Blocal{1}*Blocal{2}*Blocal{3}*Blocal{4}});


%% Prova
q00 = zeros(5,1);
pos_0 = full(B_sym(q00))*[0;0;0;1];
pos_end0 = full(B3_sym(q00))*[0;0;-l;1];
% Initialize q
Nteta = 400;
qnum = zeros(5,Nteta);
qnum(1,:) = linspace(0, 2*pi, Nteta);
% qnum(2,:) = 0.2*sin(linspace(0, 16*pi, Nteta));
qnum(3,:) = 0.7*cos(linspace(0, 16*pi, Nteta));
qnum(4,:) = 0*0.5*cos(linspace(0, 8*pi, Nteta));
qnum(5,:) = 0.1*linspace(0,Nteta,Nteta);
% Evaluate position of disk center
pos = full(B_sym(qnum));
B_d = pos;
pos_x = pos(1,4:4:end);
pos_y = pos(2,4:4:end);
pos_z = pos(3,4:4:end);
% Evaluate position of ending point pendulum
B_e = full(B3_sym(qnum));
pos_e = nan(4,Nteta);
for i = 1:Nteta
    pos_e(:,i) = B_e(:,1 + (i-1)*4: i*4)*[0;0;-l;1];
end
pos_ex = pos_e(1,:);
pos_ey = pos_e(2,:);
pos_ez = pos_e(3,:);
B_i = full(B4_sym(qnum));

figure(1)
plot3(pos_x, pos_y, pos_z)
hold on
plot3(pos_ex, pos_ey, pos_ez)
grid on
axis('equal')
% 
figure(2);
ax = gca;
t1 = hgtransform('parent',ax);
t2 = hgtransform('parent',ax);
t3 = hgtransform('parent',ax);

% Draw vertical pendulum in starting configuration (in B3 frame)
pv = line([0;pos_end0(1)],[0;pos_end0(2)],[0;pos_end0(3)],'parent', t1);
hold(ax,'on')
p1 = plot3(0,0,0,'marker','o','parent',t1);
p2 = plot3(pos_end0(1),pos_end0(2),pos_end0(3), 'marker', 'o','parent', t1);
traj = plot3(0,0,0,'color','r', 'parent', ax);

% Draw horizontal pendulum in starting configuration (in B4 frame)
ph = line([0;0],[0;l],[0;0],'parent', t2);
traj2 = plot3(0,0,0,'color','g', 'parent', ax);

% Draw disk
teta_d = linspace(0,2*pi,100) ;
disk = plot3(rd*cos(teta_d), zeros(100,1), rd*sin(teta_d), 'parent', t3);
disk_l = line([0;rd],[0;0],[0;0],'color', 'black', 'parent', t3);
disk_p = plot3(rd,0,0,'marker', 'o','color', 'r', 'parent', t3);

grid on
xlim([-2*l,2*l])
ylim([-2*l,2*l])
zlim([-2.5*l, 2*l])
view(ax, [-37.5,30])
axis('off')
pbaspect([2,2,2])
set(gcf,'color',[1,1,1]);
% for i = 1:Nteta
%     t1.Matrix = B_e(:,1 + (i-1)*4: i*4);
%     set(traj, 'XData', pos_ex(1:i), 'YData', pos_ey(1:i), 'ZData', pos_ez(1:i))
%     t2.Matrix = B_i(:,1 + (i-1)*4: i*4);
%     set(traj2, 'XData', pos_x(1:i), 'YData', pos_y(1:i), 'ZData', pos_z(1:i))
%     t3.Matrix = B_d(:,1 + (i-1)*4: i*4);
%     pause(0.02)
% end
% grid on
ph.Visible = 'off';
pv.Visible = 'off';
p1.Visible = 'off';
p2.Visible = 'off';
disk.Visible = 'off';
disk_l.Visible = 'off';
disk_p.Visible = 'off';
axis('equal')
cc = linspace(0.2,0.5,Nteta);
fill3(pos_x,pos_y,pos_z,cc)
cc2 = linspace(0.8,0.9,Nteta);
fill3(pos_ex,pos_ey,pos_ez,cc2)


%% Velocità
baseAcc = [0;0;9.81;0;0;0]; % simulating an upward gravitational acceleration

qd = SX.sym('qd', N, 1);

% twist forward propagation
twist = nan(6, N, class(q));  % twist init
twistD = nan(6, N, class(q)); % twist derivative initialization
AdGinv = cell(1, 6);          % adjoint matrices init
a = cell(6, 1);               % init relative twist ad

% first joint (track curvilinear joint) is a special case
twist(:, 1) = X(:, 1).*qd(1); % twist of the track frame
AdGinv{1} = adjointInv(Blocal{1});

for i = 2:5
    AdGinv{i} = adjointInv(Blocal{i});

    twist(:, i) = AdGinv{i}*twist(:, i-1) + X(:, i).*qd(i); % twist of i-th body (i = 1, is the base)
    a{i} = ad(X(:, i).*qd(i));
end

%% Dynamics 

tau = SX.sym('tau',5,1);

% inertia matrices
if isa(q, 'casadi.SX'); zero = casadi.SX(6,6); else; zero = zeros(6,6); end

II{1} = adjointInv(Ttz(-l/2))'*diag([mb1,mb1,mb1,(1/12)*mb1*l^2,(1/12)*mb1*l^2,0.005])*adjointInv(Ttz(-l/2));  % prima asta
II{2} = zero;  % intermedio
II{3} = zero;  % intermedio
II{4} = adjointInv(Tty(l/2))'*diag([mb2,mb2,mb2,(1/12)*mb2*l^2,0.005,(1/12)*mb2*l^2])*adjointInv(Tty(l/2)); % seconda asta
II{5} = diag([md,md,md,(1/4)*md*rd^2,(1/2)*md*rd^2,(1/4)*md*rd^2]); % disco

Mhat = cell(6, 1);             % init articulated inertia
Mbar = cell(6, 1);             % init projected inertia
bhat = nan(6, 5, class(q));    % init articulated bias
Mscalar = nan(6, 1, class(q)); % init scalarazied non structural inertia (expression collection)

for i = N:-1:1
    if i == N % final body has no children
        Mhat{i} = II{i};
        Mscalar(i) = X(:, i).'*Mhat{i}*X(:, i) + 1e-3;
        Mbar{i} = (  eye(6) - (Mhat{i}*X(:, i)*X(:, i).')./(Mscalar(i))  )*Mhat{i};

        bhat(:, i) = adStar(twist(:, i))*II{i}*twist(:, i);
        bhat(:, i) = bhat(:, i); % - Fext
        continue
    end

    j = i+1; % children body
    AdGstar = AdGinv{j}.';

    % current articulated inertia
    Mhat{i} = II{i} + AdGstar*Mbar{j}*AdGinv{j};
    Mscalar(i) = X(:, i).'*Mhat{i}*X(:, i) + 1e-3;

    % projected inertia (not required for the 1st body)
    if i ~= 1
        Mbar{i} = (  eye(6) - (Mhat{i}*X(:, i)*X(:, i).')./(Mscalar(i))  )*Mhat{i};
    end

    % current body bias
    b_i = adStar(twist(:, i))*II{i}*twist(:, i);

    % current articulation bias
    V_ij = AdGinv{j}*twist(:, i); % current body twist in children frame
    bhat(:, i) = b_i + AdGstar*...
        (...
        bhat(:, j) - Mbar{j}*a{j}*V_ij + Mhat{j}*X(:, j).*(tau(j) - X(:, j).'*bhat(:, j))./Mscalar(j)...
        ); % current articulated bias

end

Mass = Function('Mass', {[q;qd]}, {Mhat{1},Mhat{2}, Mhat{3}, Mhat{4}, Mhat{5}});

% forward propagation of accelerations
qdd = nan(5, 1, class(q));
for i = 1:N
    if i == 1 % track joint transform
        Vdot_parent = AdGinv{i}*baseAcc;
        qdd(i) = -X(:, i).'*(Mhat{i}*(Vdot_parent) + bhat(:, i))./(Mscalar(i));
        twistD(:, i) = Vdot_parent + X(:, i).*qdd(i); % base velocity is always zero
        continue
    end
    % current joint acceleration
    Vdot_parent = AdGinv{i}*twistD(:, i-1);
    V_p = a{i}*AdGinv{i}*twist(:, i-1);
    qdd(i) = (tau(i) - X(:, i).'*(Mhat{i}*(Vdot_parent - V_p) + bhat(:, i)))./(Mscalar(i));

    % current body twist derivative
    twistD(:, i) = Vdot_parent + X(:, i).*qdd(i) - V_p;

end

twistDn = Function('twistDn', {[q;qd]}, {twistD});
xdot = Function('xdot', {[q; qd], tau}, {vertcat(qd, qdd)});

%% Simulation
close all
dt = 0.01;
Nstep = 1000;
qdn = zeros(10,Nstep + 1);
qdn(end,1) = 500;%100;
c = 0*10;
for i = 1:Nstep
    u = [-0*qdn(6,i),-c*qdn(7,i),-c*qdn(8,i),-0.1*c*qdn(9,i),0]; 
    if i > 10 
        u(end) = 1;%1;
    end
    snext = RK4_step(qdn(:,i), u, xdot, dt);
    snext = full(snext);
    qdn(:,i+1) = snext;
end



% Evaluate position of disk center
q00 = zeros(5,1);
pos_0 = full(B_sym(q00))*[0;0;0;1];
pos_end0 = full(B3_sym(q00))*[0;0;-l;1];
% Evaluate position of disk center
pos = full(B_sym(qdn(1:5,:)));
B_d = pos;
pos_x = pos(1,4:4:end);
pos_y = pos(2,4:4:end);
pos_z = pos(3,4:4:end);
% Evaluate position of ending point pendulum
B_e = full(B3_sym(qdn(1:5,:)));
pos_e = nan(4,Nteta);
for i = 1:Nstep
    pos_e(:,i) = B_e(:,1 + (i-1)*4: i*4)*[0;0;-l;1];
end
pos_ex = pos_e(1,:);
pos_ey = pos_e(2,:);
pos_ez = pos_e(3,:);
B_i = full(B4_sym(qdn(1:5,:)));

figure(3)
plot3(pos_x, pos_y, pos_z)
hold on
plot3(pos_ex, pos_ey, pos_ez)
grid on
axis('equal')
% 
figure(4);
axis('equal')
ax = gca;
t1 = hgtransform('parent',ax);
t2 = hgtransform('parent',ax);
t3 = hgtransform('parent',ax);

% Draw vertical pendulum in starting configuration (in B3 frame)
line([0;pos_end0(1)],[0;pos_end0(2)],[0;pos_end0(3)],'parent', t1);
hold(ax,'on')
plot3(0,0,0,'marker','o','parent',t1);
plot3(pos_end0(1),pos_end0(2),pos_end0(3), 'marker', 'o','parent', t1);
traj = plot3(0,0,0,'color','r', 'parent', ax);

% Draw horizontal pendulum in starting configuration (in B4 frame)
line([0;0],[0;l],[0;0],'parent', t2);
traj2 = plot3(0,0,0,'color','g', 'parent', ax);

% Draw disk
teta_d = linspace(0,2*pi,100) ;
disk = plot3(rd*cos(teta_d), zeros(100,1), rd*sin(teta_d), 'parent', t3);
disk_line = line([0;rd],[0;0],[0;0],'color', 'black', 'parent', t3);
disk_point = plot3(rd,0,0,'marker', 'o','color', 'r', 'parent', t3);

grid on
xlim([-2*l,2*l])
ylim([-2*l,2*l])
zlim([-2*l-rd, 1.1*l])
view(3)
for i = 1:Nstep
    t1.Matrix = B_e(:,1 + (i-1)*4: i*4);
    set(traj, 'XData', pos_ex(1:i), 'YData', pos_ey(1:i), 'ZData', pos_ez(1:i))
    t2.Matrix = B_i(:,1 + (i-1)*4: i*4);
    set(traj2, 'XData', pos_x(1:i), 'YData', pos_y(1:i), 'ZData', pos_z(1:i))
    t3.Matrix = B_d(:,1 + (i-1)*4: i*4);
    pause(0.02)

end
grid on
