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

SelectedAxis = 4;

% Get controllers
switch SelectedAxis
    case 1
        ref_x = 3;

        mpc_x = MpcControl_x(sys_x, Ts, H);
        sys_x_start = [w(2), phi(2), v(1) , p(1)];
        %open loop
        [u_current, T_opt, X_opt, U_opt] = mpc_x.get_u(sys_x_start', ref_x);
        U_opt(:,end+1) = nan;
        plot = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys_x, xs, us); % Plot as usual
        %Close loop
        [T, X_sub, U_sub] = rocket.simulate_f( sys_x, sys_x_start, Tf, @mpc_x.get_u, ref_x);
        X = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us);
    case 2
        ref_y = 4;

        mpc_y = MpcControl_y(sys_y, Ts, H);
        sys_y_start = [w(1), phi(1), v(2) , p(2)];
        %open loop
        [u_current, T_opt, X_opt, U_opt] = mpc_y.get_u(sys_y_start', ref_y);
        U_opt(:,end+1) = nan;
        plot = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys_y, xs, us); % Plot as usual
        %Close loop
        [T, X_sub, U_sub] = rocket.simulate_f( sys_y, sys_y_start, Tf, @mpc_y.get_u, ref_y );
        X = rocket.plotvis_sub(T, X_sub, U_sub, sys_y, xs, us);
    case 3
        ref_z = 5;

        mpc_z = MpcControl_z(sys_z, Ts, H);
        sys_z_start = [v(3) , p(3)];
        %open loop
        [u_current, T_opt, X_opt, U_opt] = mpc_z.get_u(sys_z_start', ref_z);
        U_opt(:,end+1) = nan;
        plot = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys_z, xs, us); % Plot as usual
        %close loop
        [T, X_sub, U_sub] = rocket.simulate_f( sys_z, sys_z_start, Tf, @mpc_z.get_u, ref_z );
        Z = rocket.plotvis_sub(T, X_sub, U_sub, sys_z, xs, us);
    case 4
        ref_roll = deg2rad(20);

        mpc_roll = MpcControl_roll(sys_roll, Ts, H);        
        sys_roll_start = [w(3) , phi(3)];
        %open loop
        [u_current, T_opt, X_opt, U_opt] = mpc_roll.get_u(sys_roll_start', ref_roll);
        U_opt(:,end+1) = nan;
        plot = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys_roll, xs, us); % Plot as usual
        %close loop
        [T, X_sub, U_sub] = rocket.simulate_f( sys_roll, sys_roll_start, Tf, @mpc_roll.get_u, ref_roll);
        R = rocket.plotvis_sub(T, X_sub, U_sub, sys_roll, xs, us);
end

