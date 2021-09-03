classdef dms < handle
    properties
        N
        rk_steps
        ocp
        J
        eq = {}
        ieq = {}
        t0
        tf
        t
        x
        u
        p
        t0_ub
        t0_lb
        tf_ub
        tf_lb
        x_ub
        x_lb
        u_ub
        u_lb
        p_ub
        p_lb
        h
        integrals = yop.nlp_int.empty(1,0);
        timepoints = yop.nlp_timepoint.empty(1,0);
    end
    methods
        
        function obj = dms(ocp, N, rk_steps)
            obj.ocp = ocp;
            obj.N = N;
            obj.rk_steps = rk_steps;
            obj.build();
        end
        
        function obj = build(obj)
            obj.build_vars();
            obj.set_box_bnd();
            
            obj.find_timepoints_n_ints(obj.ocp.objective);
            for k=1:length(obj.ocp.eq)
                obj.find_timepoints_n_ints(obj.ocp.eq{k}.lhs);
            end
            for k=1:length(obj.ocp.ieq)
                obj.find_timepoints_n_ints(obj.ocp.ieq{k}.lhs);
            end
            
            obj.parameterize_timepoints();
            obj.parameterize_ints();
            obj.parameterize_cost();
            obj.discretize_ode();
            obj.parameterize_constraints();
        end
        
        function obj = build_vars(obj)
            % BUILD_VARS - Build nlp variable vector
            
            % Time horizon
            obj.t0 = casadi.MX.sym('t0'); %obj.ocp.t0;
            obj.tf = casadi.MX.sym('tf'); %obj.ocp.tf;
            obj.h = (obj.tf - obj.t0)/obj.N;
            
            % grid points
            obj.t = cell(obj.N+1, 1);
            obj.t{1} = obj.t0;
            obj.t{end} = obj.tf;
            for k=2:obj.N
                obj.t{k} = obj.t{k-1} + obj.h;
            end
            
            % State at grid points
            obj.x = cell(obj.N+1,1);
            for k=1:(obj.N+1)
                obj.x{k} = casadi.MX.sym(['x_' num2str(k)], obj.ocp.nx);
            end
            
            % Control intervals
            obj.u = cell(obj.N, 1);
            for k=1:obj.N
                obj.u{k} = casadi.MX.sym(['u_' num2str(k)], obj.ocp.nu);
            end
            
            % Parameters
            obj.p = casadi.MX.sym('p', obj.ocp.np);
        end
        
        function obj = set_box_bnd(obj)
            obj.t0_lb = obj.ocp.t0_lb;
            obj.t0_ub = obj.ocp.t0_ub;
            obj.tf_lb = obj.ocp.tf_lb;
            obj.tf_ub = obj.ocp.tf_ub;
            
            obj.x_lb = repmat(obj.ocp.x_lb, obj.N+1, 1);
            obj.x_ub = repmat(obj.ocp.x_ub, obj.N+1, 1);
            obj.x_lb(1:obj.ocp.nx) = obj.ocp.x0_lb;
            obj.x_ub(1:obj.ocp.nx) = obj.ocp.x0_ub;
            obj.x_lb(end-obj.ocp.nx+1:end) = obj.ocp.xf_lb;
            obj.x_ub(end-obj.ocp.nx+1:end) = obj.ocp.xf_ub;
            
            obj.u_lb = repmat(obj.ocp.u_lb, obj.N, 1);
            obj.u_ub = repmat(obj.ocp.u_ub, obj.N, 1);
            obj.u_lb(1:obj.ocp.nu) = obj.ocp.u0_lb;
            obj.u_ub(1:obj.ocp.nu) = obj.ocp.u0_ub;
            obj.u_lb(end-obj.ocp.nu+1:end) = obj.ocp.uf_lb;
            obj.u_ub(end-obj.ocp.nu+1:end) = obj.ocp.uf_ub;
                        
            obj.p_lb = obj.ocp.p_lb;
            obj.p_ub = obj.ocp.p_ub;
        end
        
        function obj = parameterize_cost(obj) 
            obj.J = obj.parameterize_expr(obj.ocp.objective);
        end
        
        function obj = discretize_ode(obj)
            F = obj.rk4();
            % Integrate trajectory individually on segments and introduce
            % continuity constraints
            for k=1:(obj.N)
                obj.add_eq( obj.x{k+1} - ...
                    F(obj.t{k}, obj.t{k+1}, obj.x{k}, obj.u{k}, obj.p) );
            end
        end
        
        function obj = parameterize_constraints(obj)
            for k=1:length(obj.ocp.eq)
                obj.add_eq(obj.parameterize_expr(obj.ocp.eq{k}.lhs));
            end
            for k=1:length(obj.ocp.ieq)
                obj.add_ieq(obj.parameterize_expr(obj.ocp.ieq{k}.lhs));
            end
        end
        
        
        function obj = find_timepoints_n_ints(obj, expr)
            [tsort, K] = topological_sort(expr);
            tsort = tsort(1:K);
            
            for k=1:K
                if isa(tsort{k}, 'yop.ast_int')
                    obj.integrals(end+1) = yop.nlp_int(tsort{k});
                elseif isa(tsort{k}, 'yop.ast_timepoint')
                    obj.timepoints(end+1) = yop.nlp_timepoint(tsort{k});
                end
            end            
            
        end
        
        function value = propagate(obj, expr)
            [tsort, K] = topological_sort(expr);
            for k=1:K
                if ~isa(tsort{k}, 'yop.ast_int') && ...
                        ~isa(tsort{k}, 'yop.ast_timepoint')
                    forward(tsort{k});
                end
            end
            % Even if last element is a timepoint or integral, it is ok to
            % evaluate, since we want to compute its value.
            value = forward(tsort{k}); 
        end
        
        
        function obj = parameterize_timepoints(obj)
            for k=1:length(obj.timepoints)
                fn = obj.ocp.expr_fn(obj.timepoints(k).expr);
                args = cell(6,1);
                [args{:}] = obj.vars_at(obj.timepoints(k).tp);
                obj.timepoints(k).value = fn(args{:});
            end
        end
        
        function obj = parameterize_ints(obj)
            % Process integrals from first to last since they are stored in
            % topological order
            
            % Compute integrand expression
            obj.ocp.set_mx();
            for k=1:length(obj.integrals)
                obj.timepoints.parameterize();
                obj.integrals.parameterize();
                value = obj.propagate(obj.integrals(k).node);
                obj.integrals(k).integrand_expr = value;
            end
            
            % Compute integrand function
            % Because it might contain timepoints, the full nlp variable
            % vector is passed as an argument.
            for k=1:length(obj.integrals)
                args = {obj.ocp.t, obj.ocp.x, obj.ocp.u, obj.ocp.p, obj.w};
                outs = {obj.integrals(k).integrand_expr};
                f = casadi.Function('f', args, outs);
                obj.integrals(k).integrand_fn = f;
            end
            
            % Integrate
            for k=1:length(obj.integrals)
                rk = obj.rk4q(obj.integrals(k).integrand_fn);
                I = 0;
                for n=1:obj.N
                    I = I + rk(obj.t{n}, obj.t{n+1}, obj.x{n}, ...
                        obj.u{n}, obj.p, obj.w);
                end
                obj.integrals(k).value = I;
            end
            
        end
        
        function val = parameterize_expr(obj, expr)
            obj.ocp.set_mx();
            if all(is_transcription_invariant(expr))
                % if the expression is transcription invariant, it means
                % that all time varying variables are dominated by a
                % timepoint, integral, or is not part of the expression.
                % This means that we can evalutate the expression, and
                % directly obtain the parametrization.
                obj.timepoints.parameterize();
                obj.integrals.parameterize();
                [tsort, K] = topological_sort(expr);
                for k=1:K
                    if ~isa(tsort{k}, 'yop.ast_int') && ...
                            ~isa(tsort{k}, 'yop.ast_timepoint')
                        forward(tsort{k});
                    end
                end
                val = value(tsort{K});
                
            else
                % If it is not invariant, it is necessary to evaluate the
                % expression for every timepoint of the problem. This is
                % done by first creating an expression function for the
                % expression, and evaluating it for every grid point.
                
                obj.timepoints.parameterize();
                obj.integrals.parameterize();
                [tsort, K] = topological_sort(expr);
                for k=1:K
                    if ~isa(tsort{k}, 'yop.ast_int') && ...
                            ~isa(tsort{k}, 'yop.ast_timepoint')
                        forward(tsort{k});
                    end
                end
                sym_expr = value(tsort{K});
                
                % Create expression function
                [tt0, ttf, tt, xx, uu, pp] = obj.ocp.get_paramlist(obj);
                args = {tt0, ttf, tt, xx, uu, pp, obj.w};
                fn = casadi.Function('int', args, {sym_expr});
                
                grid_val = cell(obj.N+1, 1);
                for n=1:obj.N
                    grid_val{n} = fn(obj.t0, obj.tf, obj.t{n}, obj.x{n}, ...
                        obj.u{n}, obj.p, obj.w);                    
                end
                grid_val{obj.N+1} = fn(obj.t0, obj.tf, obj.t{obj.N+1}, ...
                    obj.x{obj.N+1}, obj.u{obj.N}, obj.p, obj.w);
                val = vertcat(grid_val{:});
            end 
        end
        
        function rk = rk4q(obj, qfn)            
            nx = obj.ocp.nx; nu = obj.ocp.nu; np = obj.ocp.np;
            nw = numel(obj.w);
            
            % Integrator parameters
            I0 = casadi.MX.sym('i0'); % Start of integration
            If = casadi.MX.sym('if'); % End of integration
            x0 = casadi.MX.sym('x0', nx); % Initial value
            
            % ode parameters
            t_ = casadi.MX.sym('t');  % Independent variable
            x_ = casadi.MX.sym('x', nx); % State variable
            u_ = casadi.MX.sym('u', nu); % Control input
            p_ = casadi.MX.sym('p', np); % Free parameters
            w  = casadi.MX.sym('w', nw); % NLP variable vector (for tps)
            
            args = {t_,x_,u_,p_,w};
            outs = {obj.ocp.ode(t_,x_,u_,p_), qfn(t_,x_,u_,p_,w)};
            f = casadi.Function('f', args, outs);
           
            % Expression for RK4 integrator
            tt = I0;
            xx = x0;
            uu = u_;
            pp = p_;
            hh = (If-I0)/obj.rk_steps; 
            qq = 0;
            for k=1:obj.rk_steps
                [k1, q1] = f(tt     , xx        , uu, pp, w);
                [k2, q2] = f(tt+hh/2, xx+hh/2*k1, uu, pp, w);
                [k3, q3] = f(tt+hh/2, xx+hh/2*k2, uu, pp, w);
                [k4, q4] = f(tt+hh  , xx+hh  *k3, uu, pp, w);
                tt = tt + hh;
                xx = xx + hh/6*(k1 + 2*k2 + 2*k3 + k4);
                qq = qq + hh/6*(q1 + 2*q2 + 2*q3 + q4);
            end
            rk = casadi.Function('Q', {I0,If,x0,uu,pp,w}, {qq});
        end
        
        function rk = rk4(obj)            
            nx = obj.ocp.nx; nu = obj.ocp.nu; np = obj.ocp.np;
            
            % Integrator parameters
            I0 = casadi.MX.sym('I0'); % Start of integration
            If = casadi.MX.sym('If'); % End of integration
            x0 = casadi.MX.sym('x0', nx); % Initial value
            
            % ode parameters
            % (Akward naming, but otherwise Matlab's linter kills my eyes)
            t_  = casadi.MX.sym('t');  % Independent variable
            x_  = casadi.MX.sym('x', nx); % State variable
            u_  = casadi.MX.sym('u', nu); % Control input
            p_  = casadi.MX.sym('p', np); % Free parameters

            % ode rhs
            args = {t_,x_,u_,p_};
            outs = {obj.ocp.ode(t_,x_,u_,p_)};
            f = casadi.Function('f', args, outs);
           
            % Expression for RK4 integrator
            tt = I0;
            xx = x0;
            uu = u_;
            pp = p_;
            hh = (If-I0)/obj.rk_steps;
            for k=1:obj.rk_steps
                k1 = f(tt     , xx        , uu, pp);
                k2 = f(tt+hh/2, xx+hh/2*k1, uu, pp);
                k3 = f(tt+hh/2, xx+hh/2*k2, uu, pp);
                k4 = f(tt+hh  , xx+hh  *k3, uu, pp);
                tt = tt + hh;
                xx = xx + hh/6*(k1 + 2*k2 + 2*k3 + k4);
            end
            rk = casadi.Function('rk4', {I0, If, x0, uu, pp}, {xx});
        end
        
        function [t0, tf, t, x, u, p] = vars_at(obj, tp)
            t0 = obj.t0;
            tf = obj.tf;
                
            if isa(tp, 'yop.independent_initial')
                t = obj.t0;
                x = obj.x{1};
                u = obj.u{1};
                p = obj.p;
                return;
            elseif isa(tp, 'yop.ast_independent_final')
                t = obj.tf;
                x = obj.x{end};
                u = obj.u{end};
                p = obj.p;
                return;
            end
            
            [fixed, T] = fixed_horizon(obj.ocp);            
            
            if ~fixed
                error('[Yop] Problem does not have a fixed horizon')
            end
            
            dt = T/obj.N; % grid step_length    
            if yop.EQ(rem(tp,dt), 0, 1e-3)
                % grid point
                n = round(tp/dt);
                t = obj.t{n};
                x = obj.x{n};
                u = obj.u{min(obj.N, n)};
                p = obj.p;
                
            else
                % With the first statement and this floor operation n!=N+1,
                % unless t is outside the horizon, in which case an error
                % is desirable.
                % Integrate up to timepoint
                n = floor(tp/dt);
                t = tp;
                ode_int = obj.rk4();
                x = ode_int(obj.t{n}, tp, obj.x{n}, obj.u{n}, obj.p);
                u = obj.u{n};
                p = obj.p;
            end
        end
        
        function v = w(obj)
            v = vertcat(obj.t0, obj.tf, obj.x{:}, obj.u{:}, obj.p);
        end
        
        function bd = w_lb(obj)
            bd = vertcat(obj.t0_lb, obj.tf_lb, obj.x_lb, obj.u_lb, obj.p_lb);
        end
        
        function bd = w_ub(obj)
            bd = vertcat(obj.t0_ub, obj.tf_ub, obj.x_ub, obj.u_ub, obj.p_ub);
        end
        
        function n = nw(obj)
            n = length(obj.w);
        end
        
        function obj = add_eq(obj, expr)
            obj.eq{end+1} = expr;
        end
        
        function obj = add_ieq(obj, expr)
            obj.ieq{end+1} = expr;
        end
        
    end
end


% function val = disc_expr(obj, expr)
% % It is possible to have a timepoint in an integral, but not to
% % evaluate an integral at at timepoint. If the latter is
% % desired, it is necessary to use an integration state, and
% % evaluate that state at a certain timepoint. Note that Yop
% % does not test this, so evaluating an integral at a timepoint
% % is undefined behavior.
% [tsort, K] = topological_sort(expr);
% tps = yop.nlp_timepoint.empty(1,0);
% ints = yop.nlp_int.empty(1,0);
% for k=1:K
%     if isa(tsort{k}, 'yop.ast_timepoint')
%         tps(end+1) = yop.nlp_timepoint(tsort{k});
%         
%     elseif isa(tsort{k}, 'yop.ast_int')
%         ints(end+1) = yop.nlp_int(tsort{k});
%     end
% end
% 
% % Compute timepoint value.
% for k=1:length(tps)
%     fn = obj.ocp.expr_fn(tps(k).expr);
%     [tt0, ttf, ttp, xtp, utp, ptp] = obj.vars_at(tps(k).tp);
%     tps(k).value = fn(tt0, ttf, ttp, xtp, utp, ptp);
% end
% 
% % Compute integral value
% for k=1:length(ints)
%     % 1) Set the symbolic variables of the OCP to the symbolic
%     % ones, so that we can derive an expression function.
%     obj.ocp.set_sym();
%     
%     % 2) For the timepoints and integrals computed so far, they
%     % are parameterized with the computed value. If they
%     % haven't got  a value, its empty, but it should not matter
%     % as we are processing nodes in topological order, so all
%     % nodes that this node depends on has already been
%     % processed.
%     tps.parameterize();
%     ints.parameterize(); % int(expr(int(expr))), even possible??
%     
%     % 3) Forward evaluate the expression, in order to derive a
%     % parameterized expression. Integrals and timepoints are
%     % not evaluated since, we do not wish to changes their
%     % parametrization.
%     [ts, R] = topological_sort(ints(k).node);
%     for r=1:R
%         if ~isa(ts{r}, 'yop.ast_timepoint') && ...
%                 ~isa(ts{r}, 'yop.ast_int')
%             forward(ts{r});
%         end
%     end
%     integrand = forward(ts{R}); % Last is an integral, need fw.
%     
%     % 4) Now that we have an expression for the integrand, we
%     % create a casadi function from it. Beacuse it might
%     % containt timepoints, it is necessary to use the entire
%     % nlp vector as a parameter. This makes it possible to
%     % changes value of non-timepoints, while keeping the
%     % timepoint values constant.
%     
%     args = {obj.ocp.t0, obj.ocp.tf, obj.ocp.t, obj.ocp.x, ...
%         obj.ocp.u, obj.ocp.p, obj.w};
%     fn = casadi.Function('int', args, {integrand});
%     
%     % 5) Once a function for the integrand is obtained, we
%     % create a RK4 integrator for it. The integrator has the
%     % nlp vector as final parameter
%     rk4 = yop.rk4_dms_expr(obj.ocp.ode, fn, obj.ocp.nx, ...
%         obj.ocp.nu, obj.ocp.np, obj.nw, obj.rk_steps);
%     
%     % 6) Integration. The expression is integrated over every
%     % segment of the grid.
%     I = 0;
%     for n=1:obj.N
%         I = I + rk4(obj.t{n}, obj.t{n+1}, obj.t0, obj.tf, ...
%             obj.x{n}, obj.u{n}, obj.p, obj.w);
%     end
%     
%     % 7) The value of the integral is stored
%     ints(k).value = I;
%     
%     % 8) The values of the ocp variables are reset to their
%     % former value.
%     % obj.ocp.reset_sym();
% end
% 
% if all(is_transcription_invariant(expr))
%     % if the expression is transcription invariant, it means
%     % that all time varying variables are dominated by a
%     % timepoint, integral, or is not part of the expression.
%     % This means that we can evalutate the expression, and
%     % directly obtain the parametrization.
%     
%     % The symbolic variables are set in order to ensure proper
%     % execution.
%     obj.ocp.set_sym();
%     
%     % The timepoints and integrals are set to their
%     % corresponding parametrizations
%     tps.parameterize();
%     ints.parameterize();
%     
%     for k=1:K
%         tk = tsort{k};
%         if ~isa(tk,'yop.ast_timepoint')&&~isa(tk,'yop.ast_int')
%             forward(tk);
%         end
%     end
%     val = {value(tsort{K})};
%     % obj.ocp.reset_sym();
% else
%     % If it is not invariant, it is necessary to evaluate the
%     % expression for every timepoint of the problem. This is
%     % done by first creating an expression function for the
%     % expression, and evaluating it for every grid point.
%     
%     % Set symbolic values in order to derive an expression
%     % function.
%     obj.ocp.set_sym();
%     
%     % Parameterize timepoint and intergrals
%     tps.parameterize();
%     ints.parameterize();
%     
%     % Evaluate the expression
%     for k=1:K
%         tk = tsort{k};
%         if ~isa(tk,'yop.ast_timepoint')&&~isa(tk,'yop.ast_int')
%             forward(tk);
%         end
%     end
%     integrand = forward(tsort{K}); % in case timepoint or int
%     
%     % Create expression function
%     args = {obj.ocp.t0, obj.ocp.tf, obj.ocp.t, obj.ocp.x, ...
%         obj.ocp.u, obj.ocp.p, obj.w};
%     fn = casadi.Function('int', args, {integrand});
%     
%     val = cell(obj.N+1, 1);
%     for n=1:obj.N
%         val{n} = fn(obj.t0, obj.tf, obj.t{n}, obj.x{n}, ...
%             obj.u{n}, obj.p, obj.w);
%     end
%     val{obj.N+1} = fn(obj.t0, obj.tf, obj.t{obj.N+1}, ...
%         obj.x{obj.N+1}, obj.u{obj.N}, obj.p, obj.w);
%     
% end
% end