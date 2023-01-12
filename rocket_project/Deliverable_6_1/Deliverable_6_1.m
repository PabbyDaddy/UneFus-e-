 %% To Do 6.1
clc
clear
close 

addpath(fullfile('..', 'src'));

Ts = 1/20;
rocket = Rocket(Ts);

H = 2.0; %[s] Horizon length
x0 = zeros(12,1);

nmpc = NmpcControl(rocket,H);

%MPC reference with default maximum roll = 15deg
ref = @(t_,x_) ref_EPFL(t_);
Tf = 30;
rocket.anim_rate = 10; %increase this to make animation faster
[T,X,U,Ref] =rocket.simulate(x0,Tf,@nmpc.get_u,ref);
ph = rocket.plotvis(T, X,U,Ref);

% %MPC reference with specified maximum roll = 50deg
% roll_max = deg2rad(50);
% ref2 = @(t_,x_) ref_EPFL(t_,roll_max);
% Tf = 30;
% [T,X,U,Ref] =rocket.simulate(x0,Tf,@nmpc.get_u,ref2);
% ph = rocket.plotvis(T, X,U,Ref);
