%addpath(fullfile('./rocket_project/templates'));

%% TODO: This file should produce all the plots for the deliverable

w = [0, 0, 0]; %angular velocity of the body axes [wx, wy, wz]
phi=[0, 0, 0]; % Euler angle [alpha, beta, sig]
v = [0, 0, 0]; % velocities [vx, vy, vz]
p = [2, 1, 0]; % position  [x, y, z]

x0 = [deg2rad(w),deg2rad(phi), v, p]'; 

u = [deg2rad([0 0]), 0, 0 ]'; 
% (d1 d2 Pavg Pdiff) Constant input
%d1 = reflection angle of servo 1 (+/- 0.26 rad)
%d2 = reflection angle of servo 2
%Pavg = Puissance averaged  [20%, 80%]
%Pdiff= [-20% +20%]

Ts = 1/20; % Sample time
Tf = 2; %simulation time

rocket = Rocket(Ts);
[xs, us] = rocket.trim();

sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% Design MPC controller
H = 8; % Horizon length in seconds
mpc_x = MpcControl_x(sys_x, Ts, H);
% Get control input

sys_x_obj = [w(2), phi(2), v(1) , p(1)];
[u, T_opt, X_opt, U_opt] = mpc_x.get_u(sys_x_obj');
U_opt(:,end+1) = nan;

ph = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys_x, xs, us); % Plot as usual

