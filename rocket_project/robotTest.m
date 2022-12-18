clear all
close all
clc


rocket = Rocket(1/20);
Tf = 2.0; % Simulation end time

x0 = [deg2rad([0 0 0, 0 0 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state

%w = [wx, wy, wz] angular velocity of the body axes
%phi=[alpha, beta, sig] Euler angle
%v = [vx, vy, vz] velocities
%p = [x, y, z] position 

u = [deg2rad([0 0]), 60, 0 ]'; % (d1 d2 Pavg Pdiff) Constant input

%d1 = reflection angle of servo 1 (+/- 0.26 rad)
%d2 = reflection angle of servo 2
%Pavg = Puissance averaged  [20%, 80%]
%Pdiff= [-20% +20%]

[T, X, U] = rocket.simulate(x0, Tf, u); % Simulate unknown, nonlinear model
rocket.anim_rate = 1.0; % Visualize at 1.0x real−time
rocket.vis(T, X, U);


[xs, us] = rocket.trim(); % Compute steady−state for which 0 = f(xs,us)
sys = rocket.linearize(xs, us); % Linearize the nonlinear model about trim point

[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);



