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
            ubx = inf(nx, 1);
            lbx = -inf(nx, 1);
            ubu = inf(nu, 1);
            lbu = -inf(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % Cost
            %% Code section used to compute the terminal cost Qf using dlqr funtion 
            [xs, us] = rocket.trim(); 
            sys = rocket.linearize(xs, us);
            sys_d = c2d(sys,rocket.Ts);
            Q = 200*eye(nx);
            R = diag([10 10 20 10]); 
            target = [0 0 0 0 0 ref_sym(4) 0 0 0 ref_sym(1) ref_sym(2) ref_sym(3)]';
            [~,Qf,~] = dlqr(sys_d.A,sys_d.B,Q,R);
            %%
            f_symbolic = @(x,u) rocket.f(x, u);
     
            %% Compute the cost of each component with their respective weights
            cost = 0.1*sum(X_sym(1,:).^2)  + ... 
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
               20*(X_sym(:,N)-target)'*Qf*(X_sym(:,N)-target);
                    
            % Equality constraints (Casadi SX), each entry == 0

            for k=1:N-1 % loop over control intervals
            %opti.subject_to(X_sym(:,k+1) == f_discrete(X_sym(:,k), U_sym(:,k)));
            next_state = RK4(X_sym(:,k), U_sym(:,k), rocket.Ts, f_symbolic);
            opti.subject_to((X_sym(:,k+1)) ==  next_state);
            opti.subject_to(M*U_sym(:,k) <= m);
            opti.subject_to(F*X_sym(:,k) <=  f);
            %   
            end
            eq_constr = [ ; ];
            
            % Inequality constraints (Casadi SX), each entry <= 0

            % inequality constraints of the type: Gz < g
            % with G = [F zeros; zeros M) and g = [f;m]

            % F*x <= f constrain on states variables 
            F = [0 0 0 1 0 0 0 0 0 0 0 0;... %alpha
                0 0 0 -1 0 0 0 0 0 0 0 0;...
                0 0 0 0 1 0 0 0 0 0 0 0;...   %beta
                0 0 0 0 -1 0 0 0 0 0 0 0];
                        
            f = [deg2rad(85);deg2rad(85);deg2rad(85);deg2rad(85)];
           
            % M*u <= m constrians on inpute variables
            M = [1 0 0 0;... % delta 1
                 -1 0 0 0;...
                 0 1 0 0;... % delat 2
                 0 -1 0 0;...
                 0 0 1 0;... % Pavg
                 0 0 -1 0;...
                 0 0 0 1;... % Pdiff
                 0 0 0 -1];
            m = [0.26;0.26;0.26;0.26;80;-50;20;20]; % do we need to do - something for Pavg
            
            DIA_left = zeros(size(M,1), size(F, 2));
            DIA_right = zeros(size(F, 1), size(M, 2));
            
            G = [F DIA_right; DIA_left M];
            g = [f; m];
            
            ineq_constr = [G ; g];

            % For box constraints on state and input, overwrite entries of
            % lbx, ubx, lbu, ubu defined above
            
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

