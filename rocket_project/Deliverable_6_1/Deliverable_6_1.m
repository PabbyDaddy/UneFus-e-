 %% To Do 6.1
clc
clear
close 

addpath(fullfile('..', 'src'));

Ts = 1/20;
rocket = Rocket(Ts);

H = 2.0; % Horizon length in seconds
x0 = zeros(12,1);

nmpc = NmpcControl(rocket,H);

ref = @(t_,x_) ref_EPFL(t_);
Tf = 30;
rocket.anim_rate = 10; %increase this to make animation faster
[T,X,U,Ref] =rocket.simulate(x0,Tf,@nmpc.get_u,ref);
ph = rocket.plotvis(T, X,U,Ref);

