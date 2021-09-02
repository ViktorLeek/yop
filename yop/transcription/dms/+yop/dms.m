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
            obj.disc_cost();
            obj.disc_ode();
            obj.disc_con();
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
            [t0lb, t0ub, tflb, tfub] = obj.ocp.t_bd();
            obj.t0_lb = t0lb;
            obj.t0_ub = t0ub;
            obj.tf_lb = tflb;
            obj.tf_ub = tfub;
            
            [x0lb, x0ub, xlb, xub, xflb, xfub] = obj.ocp.x_bd();
            obj.x_lb = repmat(xlb, obj.N+1, 1);
            obj.x_ub = repmat(xub, obj.N+1, 1);
            obj.x_lb(1:obj.ocp.nx) = x0lb;
            obj.x_ub(1:obj.ocp.nx) = x0ub;
            obj.x_lb(end-obj.ocp.nx+1:end) = xflb;
            obj.x_ub(end-obj.ocp.nx+1:end) = xfub;
            
            [u0lb, u0ub, ulb, uub, uflb, ufub] = obj.ocp.u_bd();
            obj.u_lb = repmat(ulb, obj.N, 1);
            obj.u_ub = repmat(uub, obj.N, 1);
            obj.u_lb(1:obj.ocp.nu) = u0lb;
            obj.u_ub(1:obj.ocp.nu) = u0ub;
            obj.u_lb(end-obj.ocp.nu+1:end) = uflb;
            obj.u_ub(end-obj.ocp.nu+1:end) = ufub;
            
            [plb, pub] = obj.ocp.p_bd();
            obj.p_lb = plb;
            obj.p_ub = pub;
        end
        
        function obj = disc_cost(obj) 
            J_expr = obj.disc_expr(obj.ocp.objective);
            obj.J = J_expr{1};
            assert(length(J_expr)==1, '[Yop] Unexpected error');
        end
        
        function obj = disc_ode(obj)
            [ode_expr, ode_var] = obj.ocp.ode();
            
            F = yop.rk4_integrator(ode_expr, obj.ocp.nx, obj.ocp.nu, ...
                obj.ocp.np, obj.rk_steps);
            
            % Integrate trajectory individually on segments
            xf = cell(obj.N, 1);
            for k=1:obj.N
                xf{k} = F(obj.t{k}, obj.t{k+1}, obj.t0, obj.tf, ...
                    obj.x{k}, obj.u{k}, obj.p);
            end
            
            % Continuity
            for k=1:(obj.N-1)
                % Need to use ode_var function, in order to apply the
                % constraint on the correct states.
                x_k1 = ode_var(obj.t0, obj.tf, obj.t{k+1}, obj.x{k+1}, ...
                    obj.u{k+1}, obj.p);
                obj.add_eq( x_k1 - xf{k} );
            end
            % Last point needs special treatment as there are only N
            % control intervals
            x_k1 = ode_var(obj.t0, obj.tf, obj.t{obj.N+1}, ...
                obj.x{obj.N+1}, obj.u{obj.N}, obj.p);
            obj.add_eq( x_k1 - xf{obj.N} );
            
        end
        
        function obj = disc_con(obj)
            for k=1:length(obj.ocp.eq)
                dc = obj.disc_expr(obj.ocp.eq{k}.lhs);
                for n=1:length(dc)
                    obj.add_eq(dc{n});
                end
            end
            
            for k=1:length(obj.ocp.ieq)
                dc = obj.disc_expr(obj.ocp.ieq{k}.lhs);
                for n=1:length(dc)
                    obj.add_ieq(dc{n});
                end
            end
            
        end
        
        function val = disc_expr(obj, expr)
            % It is possible to have a timepoint in an integral, but not to
            % evaluate an integral at at timepoint. If the latter is
            % desired, it is necessary to use an integration state, and
            % evaluate that state at a certain timepoint. Note that Yop
            % does not test this, so evaluating an integral at a timepoint
            % is undefined behavior.
            [tsort, K] = topological_sort(expr);
            tps = yop.nlp_timepoint.empty(1,0); 
            ints = yop.nlp_int.empty(1,0);
            for k=1:K
                if isa(tsort{k}, 'yop.ast_timepoint')
                    tps(end+1) = yop.nlp_timepoint(tsort{k});
                    
                elseif isa(tsort{k}, 'yop.ast_int')
                    ints(end+1) = yop.nlp_int(tsort{k});
                end
            end
            
            % Compute timepoint value.
            for k=1:length(tps)
                fn = obj.ocp.expr_fn(tps(k).expr);
                [tt0, ttf, ttp, xtp, utp, ptp] = obj.vars_at(tps(k).tp);
                tps(k).value = fn(tt0, ttf, ttp, xtp, utp, ptp);
            end
            
            % Compute integral value
            for k=1:length(ints)
                % 1) Set the symbolic variables of the OCP to the symbolic
                % ones, so that we can derive an expression function.
                obj.ocp.set_sym();
                
                % 2) For the timepoints and integrals computed so far, they
                % are parameterized with the computed value. If they
                % haven't got  a value, its empty, but it should not matter
                % as we are processing nodes in topological order, so all
                % nodes that this node depends on has already been
                % processed.
                tps.parameterize();
                ints.parameterize(); % int(expr(int(expr))), even possible??
                
                % 3) Forward evaluate the expression, in order to derive a
                % parameterized expression. Integrals and timepoints are
                % not evaluated since, we do not wish to changes their
                % parametrization.
                [ts, R] = topological_sort(ints(k).node);
                for r=1:R
                    if ~isa(ts{r}, 'yop.ast_timepoint') && ...
                            ~isa(ts{r}, 'yop.ast_int')
                        forward(ts{r});
                    end
                end
                integrand = forward(ts{R}); % Last is an integral, need fw.
                
                % 4) Now that we have an expression for the integrand, we
                % create a casadi function from it. Beacuse it might
                % containt timepoints, it is necessary to use the entire
                % nlp vector as a parameter. This makes it possible to
                % changes value of non-timepoints, while keeping the
                % timepoint values constant.

                args = {obj.ocp.t0, obj.ocp.tf, obj.ocp.t, obj.ocp.x, ...
                    obj.ocp.u, obj.ocp.p, obj.w};
                fn = casadi.Function('int', args, {integrand});
                
                % 5) Once a function for the integrand is obtained, we
                % create a RK4 integrator for it. The integrator has the
                % nlp vector as final parameter
                rk4 = yop.rk4_dms_expr(obj.ocp.ode, fn, obj.ocp.nx, ...
                    obj.ocp.nu, obj.ocp.np, obj.nw, obj.rk_steps);
                
                % 6) Integration. The expression is integrated over every
                % segment of the grid.
                I = 0;
                for n=1:obj.N
                    I = I + rk4(obj.t{n}, obj.t{n+1}, obj.t0, obj.tf, ...
                        obj.x{n}, obj.u{n}, obj.p, obj.w);
                end
                
                % 7) The value of the integral is stored
                ints(k).value = I;
                
                % 8) The values of the ocp variables are reset to their
                % former value.
                % obj.ocp.reset_sym();
            end
            
            if all(is_transcription_invariant(expr))
                % if the expression is transcription invariant, it means
                % that all time varying variables are dominated by a
                % timepoint, integral, or is not part of the expression.
                % This means that we can evalutate the expression, and
                % directly obtain the parametrization.
                
                % The symbolic variables are set in order to ensure proper
                % execution. 
                obj.ocp.set_sym();
                
                % The timepoints and integrals are set to their
                % corresponding parametrizations
                tps.parameterize();
                ints.parameterize();
                
                for k=1:K
                    tk = tsort{k};
                    if ~isa(tk,'yop.ast_timepoint')&&~isa(tk,'yop.ast_int')
                        forward(tk);
                    end
                end
                val = {value(tsort{K})};
                % obj.ocp.reset_sym();
            else
                % If it is not invariant, it is necessary to evaluate the
                % expression for every timepoint of the problem. This is
                % done by first creating an expression function for the
                % expression, and evaluating it for every grid point.
                
                % Set symbolic values in order to derive an expression
                % function.
                obj.ocp.set_sym();
                
                % Parameterize timepoint and intergrals
                tps.parameterize();
                ints.parameterize();
                
                % Evaluate the expression
                for k=1:K
                    tk = tsort{k};
                    if ~isa(tk,'yop.ast_timepoint')&&~isa(tk,'yop.ast_int')
                        forward(tk);
                    end
                end
                integrand = forward(tsort{K}); % in case timepoint or int
                
                % Create expression function
                args = {obj.ocp.t0, obj.ocp.tf, obj.ocp.t, obj.ocp.x, ...
                    obj.ocp.u, obj.ocp.p, obj.w};
                fn = casadi.Function('int', args, {integrand});
                
                val = cell(obj.N+1, 1);
                for n=1:obj.N
                    val{n} = fn(obj.t0, obj.tf, obj.t{n}, obj.x{n}, ...
                        obj.u{n}, obj.p, obj.w);                    
                end
                val{obj.N+1} = fn(obj.t0, obj.tf, obj.t{obj.N+1}, ...
                    obj.x{obj.N+1}, obj.u{obj.N}, obj.p, obj.w);
                
            end            
        end
        
        function [t0, tf, t, x, u, p] = vars_at(obj, tp)
            t0 = obj.t0;
            tf = obj.tf;
                
            if isa(tp, 'yop.independent_initial')
                t = obj.t{1};
                x = obj.x{1};
                u = obj.u{1};
                p = obj.p;
                return;
            elseif isa(tp, 'yop.ast_independent_final')
                t = obj.t{end};
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
            if abs(rem(tp, dt)) <= 1e-9
                % No integration required as the point is close enough to a
                % grid point.
                n = round(tp/dt);
                t = obj.t{n};
                x = obj.x{n};
                if n == obj.N+1
                    u = obj.u{obj.N};
                else
                    u = obj.u{n};
                end
                p = obj.p;
                
            else
                % With the first statement and this floor operation n!=N+1,
                % unless t is outside the horizon, in which case an error
                % is desirable.
                t = tp;
                ode_int = yop.rk4_integrator(obj.ocp.ode, obj.ocp.nx, ...
                    obj.ocp.nu, obj.ocp.np, obj.rk_steps);
                n = floor(tp/dt);
                x = ode_int(obj.t{n}, tp, obj.t0, obj.tf, ...
                    obj.x{n}, obj.u{n}, obj.p);
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