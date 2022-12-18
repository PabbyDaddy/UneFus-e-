%addpath(fullfile('./rocket_project/templates'));

%% TODO: This file should produce all the plots for the deliverable

x0 = [deg2rad([0 0 0, 0 0 0]), 0 0 0, 2 1 0]'; 
% (w, phi, v, p) Initial state
%w = [wx, wy, wz] angular velocity of the body axes
%phi=[alpha, beta, sig] Euler angle
%v = [vx, vy, vz] velocities
%p = [x, y, z] position 

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

% x = [ωy , β, vx , x]
[u, T_opt, X_opt, U_opt] = mpc_x.get_u([0 0 0 2]');
U_opt(:,end+1) = nan;

ph = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys_x, xs, us); % Plot as usual

