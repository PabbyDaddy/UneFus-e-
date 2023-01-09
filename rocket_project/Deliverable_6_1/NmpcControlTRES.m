function opti_eval = NmpcControlTRES(rocket, H)

import casadi.*
opti = casadi.Opti(); % Optimization problem

N = ceil(H/rocket.Ts); % MPC horizon
nx = 12; % Number of states
nu = 4;  % Number of inputs

% Decision variables (symbolic)
X_sym = opti.variable(nx, N); % state trajectory
U_sym = opti.variable(nu, N-1);   % control trajectory)

% Parameters (symbolic)
x0_sym  = opti.parameter(nx, 1);  % initial state
ref_sym = opti.parameter(4, 1);   % target position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
% h = rocket.Ts; %Just for simplicity
% k1 = @(x,u) rocket.f(x,u);
% k2 = @(x,u) rocket.f(x+h/2*k1(x,u), u);
% k3 = @(x,u) rocket.f(x+h/2*k2(x,u), u);
% k4 = @(x,u) rocket.f(x+h*k3(x,u),   u);
% f_discrete = @(x,u) x + h/6*(k1(x,u)+2*k2(x,u)+2*k3(x,u)+k4(x,u));
% f_ano = @(x,u) rocket.f(x,u);
% 
% opti.subject_to(X_sym(:,1)==x0_sym);   % use initial position
% 
% for i = 1:(N-1)
%     opti.subject_to(X_sym(i+1)==RK4(X_sym(:,i), U_sym(:,i), rocket.Ts, f_ano));
% end
% 
% opti.subject_to(-deg2rad(85) <= X_sym(4,:));% Limit on alpha
% opti.subject_to(X_sym(4,:) <= deg2rad(85));% Limit on alpha
% opti.subject_to(-deg2rad(85) <= X_sym(5,:));% Limit on beta
% opti.subject_to(X_sym(5,:) <= deg2rad(85));% Limit on beta
% opti.subject_to(-deg2rad(85) <= X_sym(6,:));% Limit on gamma
% opti.subject_to(X_sym(6,:) <= deg2rad(85));% Limit on gamma
% 
% opti.subject_to(-deg2rad(15) <= U_sym(1,:)); %Limit on d1
% opti.subject_to(U_sym(1,:) <= deg2rad(15)); %Limit on d1
% opti.subject_to(-deg2rad(15) <= U_sym(2,:)); %Limit on d2
% opti.subject_to(U_sym(2,:) <= deg2rad(15)); %Limit on d2
% opti.subject_to(50 <= U_sym(3,:)); %Limit on Pavg
% opti.subject_to(U_sym(3,:) <= 80); %Limit on Pavg
% opti.subject_to(-20 <= U_sym(4,:));%Limit on Pdiff
% opti.subject_to(U_sym(4,:) <= 20);%Limit on Pdiff
% 
% 
% opti.minimize(...
%     0.1*sum(X_sym(1,:).^2) + ...
%     0.1*sum(X_sym(2,:).^2) + ...
%     0.1*sum(X_sym(3,:).^2) + ...
%     0.1*sum(X_sym(4,:).^2) + ...
%     0.1*sum(X_sym(5,:).^2) + ...
%     0.1*sum(X_sym(7,:).^2) + ...
%     0.1*sum(X_sym(8,:).^2) + ...
%     0.1*sum(X_sym(9,:).^2) + ...
%     100*sum((X_sym(10,:) - ref_sym(1,:)).^2)  + ... x ref
%     100*sum((X_sym(11,:) - ref_sym(2,:)).^2) + ... y ref
%     100*sum((X_sym(12,:) - ref_sym(3,:)).^2) + ... z ref
%     100*sum((X_sym(6,:) - ref_sym(4,:)).^2) + ... roll ref
%     U_sym(1,:)*U_sym(1,:)' + U_sym(2,:)*U_sym(2,:)' + U_sym(3,:)*U_sym(3,:)' + U_sym(4,:)*U_sym(4,:)');



[xs, us] = rocket.trim(); 
sys = rocket.linearize(xs, us);
sys_d = c2d(sys,rocket.Ts);
Q = 200*eye(nx);
R = diag([10 10 20 10]); 
target = [0 0 0 0 0 ref_sym(4) 0 0 0 ref_sym(1) ref_sym(2) ref_sym(3)]';
[~,Qf,~] = dlqr(sys_d.A,sys_d.B,Q,R);
% F*x <= f
F = [0 0 0 1 0 0 0 0 0 0 0 0;... %alpha
    0 0 0 -1 0 0 0 0 0 0 0 0;...
    0 0 0 0 1 0 0 0 0 0 0 0;...   %beta
    0 0 0 0 -1 0 0 0 0 0 0 0];
            
f = [deg2rad(85);deg2rad(85);deg2rad(85);deg2rad(85)];

% M*u <= m
M = [1 0 0 0;... % delta 1
     -1 0 0 0;...
     0 1 0 0;... % delat 2
     0 -1 0 0;...
     0 0 1 0;... % Pavg
     0 0 -1 0;...
     0 0 0 1;... % Pdiff
     0 0 0 -1];
m = [0.26;0.26;0.26;0.26;80;-50;20;20]; % do we need to do - something for Pavg

%f_discrete = @(x,u) RK4(x,u,rocket.Ts,rocket.f(x, u));
f_symbolic = @(x,u) rocket.f(x, u);
% ---- objective ---------
 opti.minimize(...
   0.1*sum(X_sym(1,:).^2)  + ... 
   0.1*sum(X_sym(2,:).^2)  + ... 
   0.1*sum(X_sym(3,:).^2)  + ... 
   0.1*sum(X_sym(4,:).^2)  + ... 
   0.1*sum(X_sym(5,:).^2)  + ... 
   200*( sum( (X_sym(6,:)-ref_sym(4)).^2 ))  + ... 
   0.1*sum(X_sym(7,:).^2)  + ... 
   0.1*sum(X_sym(8,:).^2)  + ... 
   0.1*sum(X_sym(9,:).^2)  + ... 
   200*( sum( (X_sym(10,:)-ref_sym(1)).^2 ))  + ... 
   200*( sum( (X_sym(11,:)-ref_sym(2)).^2 ))  + ... 
   200*( sum( (X_sym(12,:)-ref_sym(3)).^2 ))  + ... 
   0.05*U_sym(1,:)*U_sym(1,:)' + ... 
   0.05*U_sym(2,:)*U_sym(2,:)' + ...
   1*U_sym(3,:)*U_sym(3,:)' + ...
   0.01*U_sym(4,:)*U_sym(4,:)'+ ...
   20*(X_sym(:,N)-target)'*Qf*(X_sym(:,N)-target));


opti.subject_to(X_sym(:,1)==x0_sym);

 for k=1:N-1 % loop over control intervals
%   opti.subject_to(X_sym(:,k+1) == f_discrete(X_sym(:,k), U_sym(:,k)));
    next_state = RK4(X_sym(:,k), U_sym(:,k), rocket.Ts, f_symbolic);
    opti.subject_to((X_sym(:,k+1)) ==  next_state);
    opti.subject_to(M*U_sym(:,k) <= m);
    opti.subject_to(F*X_sym(:,k) <=  f);
%   
 end


% 
% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Setup solver ------
ops = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
opti.solver('ipopt', ops);

% Create function to solve and evaluate opti
opti_eval = @(x0_, ref_) solve(x0_, ref_, opti, x0_sym, ref_sym, U_sym);
end

function u = solve(x0, ref, opti, x0_sym, ref_sym, U_sym)

% ---- Set the initial state and reference ----
opti.set_value(x0_sym, x0);
opti.set_value(ref_sym, ref);

% ---- Solve the optimization problem ----
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U_sym(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end
