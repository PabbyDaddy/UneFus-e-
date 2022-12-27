%addpath(fullfile('./rocket_project/templates'));

%% TODO: This file should produce all the plots for the deliverable

w = [0, 0, 0]; %angular velocity of the body axes [wx, wy, wz]
phi= deg2rad([0, 0, 0]); % Euler angle [alpha, beta, sig]
v = deg2rad([0, 0, 0]); % velocities [vx, vy, vz]
p = [0, 0, 0]; % position  [x, y, z]

x0 = [w,phi, v, p]'; 

u0 = [deg2rad([0 0]), 0, 0 ]'; 
% (d1 d2 Pavg Pdiff) Constant input
%d1 = reflection angle of servo 1 (+/- 0.26 rad)
%d2 = reflection angle of servo 2
%Pavg = Puissance averaged  [20%, 80%]
%Pdiff= [-20% +20%]

Ts = 1/20; % Sample time
Tf = 3; %simulation time

rocket = Rocket(Ts);
[xs, us] = rocket.trim();

sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

% Design MPC controller
H = 8; % Horizon length in seconds



mpc_x = MpcControl_x(sys_x, Ts, H);
mpc_y = MpcControl_y(sys_y, Ts, H);
mpc_z = MpcControl_z(sys_z, Ts, H);
mpc_roll = MpcControl_roll(sys_roll, Ts, H);

mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
% Evaluate once and plot optimal open−loop trajectory,
% pad last input to get consistent size with time and state

%define the ref to go
ref4 = [2 2 2 deg2rad(40)]';

[u, T_opt, X_opt, U_opt] = mpc.get_u(x0, ref4);
U_opt(:,end+1) = nan;

ph = rocket.plotvis(T_opt, X_opt, U_opt, ref4); % Plot as usual
ph.fig.Name = 'MPC in linear simulation'; % Set a figure title



% Setup reference function
ref = @(t_ , x ) ref_EPFL(t_);
% Simulate
Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
% Visualize
rocket.anim_rate = 2; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title
