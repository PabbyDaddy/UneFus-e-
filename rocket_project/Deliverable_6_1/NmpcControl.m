classdef NmpcControl < handle
    
    properties
        solver
        nx, nu, N
        nlp_x0
        nlp_lbx, nlp_ubx
        nlp_lbg, nlp_ubg
        nlp_p
        
        T_opt
        sol
        idx
        
        % Warmstart
        nlp_lam_x0
        nlp_lam_g0
    end
    
    methods
        function obj = NmpcControl(rocket, tf)
           
            import casadi.*
            
            N_segs = ceil(tf/rocket.Ts); % MPC horizon
            nx = 12; % Number of states
            nu = 4;  % Number of inputs
            
            % Decision variables (symbolic)
            N = N_segs + 1; % Index of last point
            X_sym = SX.sym('X_sym', nx, N); % state trajectory
            U_sym = SX.sym('U_sym', nu, N-1); % control trajectory)
            
            % Parameters (symbolic)
            x0_sym  = SX.sym('x0_sym', nx, 1);  % initial state
            ref_sym = SX.sym('ref_sym', 4, 1);  % target position
            
            % Default state and input constraints
            ubx = inf(nx, 1); %Upper bound x
            lbx = -inf(nx, 1); %Lower bound x
            ubu = inf(nu, 1); %Upper bound u
            lbu = -inf(nu, 1); %Lower bound u
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            opti = casadi.Opti(); % Optimization problem
            

            %% Code section used to compute the terminal cost Qf using dlqr funtion 
            [xs, us] = rocket.trim(); %steady state valus for x and u
            sys = rocket.linearize(xs, us); %linearize the rocket around xs and us
            sys_d = c2d(sys,rocket.Ts);
            Q = 300*eye(nx); %state cost 
            R = diag([10 10 20 10]); %inpute cost
            target = [0 0 0 0 0 ref_sym(4) 0 0 0 ref_sym(1) ref_sym(2) ref_sym(3)]'; %target (EPFL logo)
            [~,Qf,~] = dlqr(sys_d.A,sys_d.B,Q,R); %compute the infinte lqr horizon cost 
            f_symbolic = @(x_, u_) rocket.f(x_,u_);

            %% Compute the cost

            %Tuning on the weights W1, W2, W3, W4, W5, W6, W7
            W1 = 20; % 10 %5 %NW 10
            W2 = 100; %50 % 100 %NW 50
            W3 = 10; %1
            W4 = 200; %100 %200 %NW 100
            W5 = 0.1; % 0.01 %0.05 %NW 0.5
            W6 = 0.01; % 0.0001 %0.005 %NW 0.05
            W7 = 2; %2 %NW 5

            offset = 0;%56.667;

           cost = W1*sum(X_sym(1,:).^2)  + ... 
               W1*sum(X_sym(2,:).^2)  + ... 
               W1*sum(X_sym(3,:).^2)  + ... 
               W1*sum(X_sym(4,:).^2)  + ... 
               W1*sum(X_sym(5,:).^2)  + ... 
               W2*sum( (X_sym(6,:)-ref_sym(4)).^2 )  + ... 
               W3*sum(X_sym(7,:).^2)  + ... 
               W3*sum(X_sym(8,:).^2)  + ... 
               W3*sum(X_sym(9,:).^2)  + ... 
               W2*sum( (X_sym(10,:)-ref_sym(1)).^2 )  + ... 
               W2*sum( (X_sym(11,:)-ref_sym(2)).^2 )  + ... 
               W4*sum( (X_sym(12,:)-ref_sym(3)).^2 )  + ... 
               W5*sum(U_sym(1,:).^2) + ... 
               W5*sum(U_sym(2,:).^2) + ...
               W6*sum((U_sym(3,:) - offset).^2) + ...
               W6*sum(U_sym(4,:).^2)+ ...
               W7*(X_sym(:,N)-target)'*Qf*(X_sym(:,N)-target);

           
            %% Inequality constraints (Casadi SX), each entry <= 0

               % F*x <= f constrain on states variables 
            F = [0 0 0 1 0 0 0 0 0 0 0 0;... % alpha upper bound
                0 0 0 -1 0 0 0 0 0 0 0 0;... % alpha lower bound
                0 0 0 0 1 0 0 0 0 0 0 0;...   % beta upper bound
                0 0 0 0 -1 0 0 0 0 0 0 0];   %beta lower bound
                        
            f = [ deg2rad(80) ; deg2rad(80) ; deg2rad(80) ; deg2rad(80) ];
           
            % M*u <= m constrians on inpute variables
            M = [1 0 0 0;... % delta 1 upper bound
                 -1 0 0 0;... % delta 1 lower bound
                 0 1 0 0;... % delta 2 upper bound
                 0 -1 0 0;...% delta 2 lower bound
                 0 0 1 0;... % Pavg upper bound
                 0 0 -1 0;...% Pavg lower bound
                 0 0 0 1;... % Pdiff upper bound
                 0 0 0 -1]; % Pdiff lower bound

            m = [ deg2rad(15) ; deg2rad(15) ; deg2rad(15) ; deg2rad(15) ; 80; -50 ; 20 ; 20 ]; 

            Gf = (F*X_sym - f);
            Gu = (M*U_sym - m);

            ineq_constr = [  Gf(1, 1:(end)), Gf(2, 1:(end)), Gf(3, 1:(end)), Gf(4, 1:(end)),...
                 Gu(1, (1:end)), Gu(2, (1:end)), Gu(3, (1:end)), Gu(4, (1:end)), Gu(5, (1:end)),...
                 Gu(6, (1:end)), Gu(7, (1:end)), Gu(8, (1:end))  ]';

             %% Equality constraints (Casadi SX), each entry == 0
            
            eq_constr = (X_sym(:,1) - x0_sym);

            for k=1:N-1 % loop over control intervals
            next_state = RK4(X_sym(:,k), U_sym(:,k), rocket.Ts, f_symbolic);
            eq_constr = [ eq_constr; ((X_sym(:,k+1)) -  next_state) ];
            end

            %% 
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ---- Assemble NLP ------
            nlp_x = [X_sym(:); U_sym(:)];
            nlp_p = [x0_sym; ref_sym];
            nlp_f = cost;
            nlp_g = [eq_constr; ineq_constr];
            
            nlp = struct('x', nlp_x, 'p', nlp_p, 'f', nlp_f, 'g', nlp_g);
            
            % ---- Setup solver ------
            opts = struct('ipopt', struct('print_level', 0), 'print_time', false);
            obj.solver = nlpsol('solver', 'ipopt', nlp, opts);
            
            % ---- Assemble NLP bounds ----
            obj.nlp_x0  = zeros(size(nlp_x));
            
            obj.nlp_ubx = [repmat(ubx, N, 1); repmat(ubu, (N-1), 1)];
            obj.nlp_lbx = [repmat(lbx, N, 1); repmat(lbu, (N-1), 1)];
            
            obj.nlp_ubg = [zeros(size(eq_constr)); zeros(size(ineq_constr))];
            obj.nlp_lbg = [zeros(size(eq_constr)); -inf(size(ineq_constr))];
            
            obj.nlp_p = [zeros(size(x0_sym)); zeros(size(ref_sym))];
            
            obj.nlp_lam_x0 = [];
            obj.nlp_lam_g0 = [];
            
            obj.nx = nx;
            obj.nu = nu;
            obj.N = N;
            obj.T_opt = linspace(0, N * rocket.Ts, N);
            
            obj.idx.X = [1, obj.N * obj.nx];
            obj.idx.U = obj.idx.X(2) + [1, (obj.N-1) * obj.nu];
            obj.idx.u0 = obj.idx.U(1) + [0, obj.nu-1];
        end
        
        function [u, T_opt, X_opt, U_opt] = get_u(obj, x0, ref)
            
            obj.solve(x0, ref);
            
            % Evaluate u0
            nlp_x = obj.sol.x;
            id = obj.idx.u0;
            u = full( nlp_x(id(1):id(2)) );      
            
            if nargout > 1, T_opt = obj.get_T_opt(); end
            if nargout > 2, X_opt = obj.get_X_opt(); end
            if nargout > 3, U_opt = obj.get_U_opt(); end
            return
            
            % Additional evaluation
            % Complete trajectory
            % % X_opt = full(reshape(nlp_x(idx_X(1):idx_X(2)), obj.nx, obj.N));
            % % U_opt = full(reshape(nlp_x(idx_U(1):idx_U(2)), obj.nu, obj.N - 1));
            % %
            % % cost_opt = full(sol.f);
            % % constr_opt = full(sol.g);
            % %
            % % stats = obj.solver.stats;
        end
        
        function solve(obj, x0, ref)
            
            % ---- Set the initial state and reference ----
            obj.nlp_p = [x0; ref];     % Initial condition
            obj.nlp_x0(1:obj.nx) = x0; % Initial guess consistent
            
            % ---- Solve the optimization problem ----
            args = {'x0', obj.nlp_x0, ...
                'lbg', obj.nlp_lbg, ...
                'ubg', obj.nlp_ubg, ...
                'lbx', obj.nlp_lbx, ...
                'ubx', obj.nlp_ubx, ...
                'p', obj.nlp_p, ...
                %                 'lam_x0', obj.nlp_lam_x0, ...
                %                 'lam_g0', obj.nlp_lam_g0
                };
            
            obj.sol = obj.solver(args{:});
            if obj.solver.stats.success ~= true
                solve_status_str = obj.solver.stats.return_status;
                fprintf([' [' class(obj) ': ' solve_status_str '] ']);
                obj.sol.x(obj.idx.u0) = nan;
            end
            
            % Use the current solution to speed up the next optimization
            obj.nlp_x0 = obj.sol.x;
            obj.nlp_lam_x0 = obj.sol.lam_x;
            obj.nlp_lam_g0 = obj.sol.lam_g;
        end
        function T_opt = get_T_opt(obj)
            T_opt = obj.T_opt;
        end
        function X_opt = get_X_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.X;
            X_opt = full(reshape(nlp_x(id(1):id(2)), obj.nx, obj.N));
        end
        function U_opt = get_U_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.U;
            U_opt = full(reshape(nlp_x(id(1):id(2)), obj.nu, obj.N - 1));
        end
    end
end
