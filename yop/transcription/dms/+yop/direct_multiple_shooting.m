classdef direct_multiple_shooting < handle
    properties
        N
        rk_steps
        t0
        tf
        t
        x
        u
        p
    end
    methods
        
        function obj = direct_multiple_shooting(N, rk_steps)
            obj.N = N;
            obj.rk_steps = rk_steps;
        end
        
        function sol = solve(obj, ocp)
            obj.build_vars(ocp);
            [w_lb, w_ub] = obj.set_box_bnd(ocp);
            
            tps = obj.parameterize_timepoints(ocp);
            ints = obj.parameterize_integrals(ocp, tps);
            J = obj.parameterize_objective(ocp, tps, ints);
            g = obj.parameterize_equality(ocp, tps, ints);
            h = obj.parameterize_inequality(ocp, tps, ints);
            cc = obj.discretize_ode(ocp); % Continuity constraints
            
            sol = obj.solve_nlp(J, [g; cc], h, w_lb, w_ub);
        end
        
        function sol = solve_nlp(obj, J, g, h, w_lb, w_ub) %, options
            w0 = ones(size(obj.w));
            nlp.f = J;
            nlp.x = obj.w;
            nlp.g = [g; h];
            solver = casadi.nlpsol('solver', 'ipopt', nlp);
            sol = solver( ...
                'x0', w0, ...
                'lbx', w_lb, ...
                'ubx', w_ub, ...
                'ubg', [zeros(size(g)); zeros(size(h))], ...
                'lbg', [zeros(size(g)); -inf(size(h))] ...
                );
        end
        
        function obj = build_vars(obj, ocp)
            % BUILD_VARS - Build nlp variable vector
            
            % Time horizon
            obj.t0 = casadi.MX.sym('t0');
            obj.tf = casadi.MX.sym('tf');
            
            % grid points
            obj.t = cell(obj.N+1, 1);
            obj.t{1} = obj.t0;
            obj.t{end} = obj.tf;
            dt = (obj.tf-obj.t0)/obj.N;
            for k=2:obj.N
                obj.t{k} = obj.t{k-1} + dt;
            end
            
            % State at grid points
            obj.x = cell(obj.N+1,1);
            for k=1:(obj.N+1)
                obj.x{k} = casadi.MX.sym(['x_' num2str(k)], ocp.n_x);
            end
            
            % Control intervals
            obj.u = cell(obj.N, 1);
            for k=1:obj.N
                obj.u{k} = casadi.MX.sym(['u_' num2str(k)], ocp.n_u);
            end
            
            % Parameters
            obj.p = casadi.MX.sym('p', ocp.n_p);
        end
        
        function [w_lb, w_ub] = set_box_bnd(obj, ocp)
            t0_lb = ocp.t0_lb;
            t0_ub = ocp.t0_ub;
            tf_lb = ocp.tf_lb;
            tf_ub = ocp.tf_ub;
            
            x_lb = repmat(ocp.x_lb, obj.N+1, 1);
            x_ub = repmat(ocp.x_ub, obj.N+1, 1);
            x_lb(1 : ocp.n_x) = ocp.x0_lb;
            x_ub(1 : ocp.n_x) = ocp.x0_ub;
            x_lb(end - ocp.n_x + 1 : end) = ocp.xf_lb;
            x_ub(end - ocp.n_x + 1 : end) = ocp.xf_ub;
            
            u_lb = repmat(ocp.u_lb, obj.N, 1);
            u_ub = repmat(ocp.u_ub, obj.N, 1);
            u_lb(1 : ocp.n_u) = ocp.u0_lb;
            u_ub(1 : ocp.n_u) = ocp.u0_ub;
            u_lb(end - ocp.n_u + 1 : end) = ocp.uf_lb;
            u_ub(end - ocp.n_u + 1 : end) = ocp.uf_ub;
                        
            p_lb = ocp.p_lb;
            p_ub = ocp.p_ub;
            
            w_ub = vertcat(t0_ub, tf_ub, x_ub, u_ub, p_ub);
            w_lb = vertcat(t0_lb, tf_lb, x_lb, u_lb, p_lb);
        end
        
        function cc = discretize_ode(obj, ocp)
            F = obj.rk4(ocp);
            % Integrate trajectory individually on segments and introduce
            % continuity constraints
            cc = [];
            for k=1:(obj.N)
                ck = obj.x{k+1} - ...
                    F(obj.t{k}, obj.t{k+1}, obj.x{k}, obj.u{k}, obj.p);
                cc = [cc(:); ck(:)];
            end
        end
        
        
        function tps = parameterize_timepoints(obj, ocp)
            % Important that this is in topological order, should be
            % handled by ocp. The reason is that they need to be
            % parameterized in order, otherwise there is no gurantee of the
            % results. As more and more timepoints are discretized they can
            % start appearing in expressions with other timepoints, and
            % since this is done in topological order, there is no
            % conflict.
            ints = zeros(ocp.n_int, 1);
            tps = [];
            for k = 1:length(ocp.timepoints)
                fn = ocp.timepoints(k).fn;
                [tt,xx,uu,pp] = ...
                    obj.vars_at(ocp.timepoints(k).ast.timepoint);
                tmp = [tps; zeros(ocp.n_tp-length(tps),1)];
                tp_k = fn(tt, xx, uu, pp, tmp, ints);
                tps = [tps(:); tp_k(:)];
            end
        end
        
        function ints = parameterize_integrals(obj, ocp, tps)
            % Important that this is in topological order, should be
            % handled by ocp. The reason is that they need to be
            % parameterized in order, otherwise there is no gurantee of the
            % results. As more and more timepoints are discretized they can
            % start appearing in expressions with other timepoints, and
            % since this is done in topological order, there is no
            % conflict.
            ints = [];
            for k = 1:length(ocp.integrals)
                rk = obj.rk4q(ocp, ocp.integrals(k).fn);
                tmp = [ints; zeros(ocp.n_int-length(ints),1)];
                I = 0;
                for n=1:obj.N
                    I = I + rk(obj.t{n}, obj.t{n+1}, obj.x{n}, ...
                        obj.u{n}, obj.p, tps, tmp);
                end
                ints = [ints(:); I(:)];
            end
        end
        
        function J = parameterize_objective(obj, ocp, tps, ints)
            J = obj.parameterize_expression(ocp.objective, tps, ints);
        end
        
        function g = parameterize_equality(obj, ocp, tps, ints)
            g = [];
            for k=1:length(ocp.equality_constraints)
                g_k = obj.parameterize_expression( ...
                    ocp.equality_constraints(k), tps, ints);
                g = [g(:); g_k(:)];
            end
        end
        
        function h = parameterize_inequality(obj, ocp, tps, ints)
            h = [];
            for k=1:length(ocp.inequality_constraints)
                h_k = obj.parameterize_expression( ...
                    ocp.inequality_constraints(k), tps, ints);
                h = [h(:); h_k(:)];
            end
        end
        
        function disc = parameterize_expression(obj, e, tps, ints)
            if all(is_transcription_invariant(e.ast))
                tmp = e.fn(obj.t{1},obj.x{1},obj.u{1},obj.p,tps, ints);
                disc = tmp(:);
                
            else
                disc = [];
                for n=1:obj.N
                    tmp = e.fn(obj.t{n},obj.x{n},obj.u{n},obj.p,tps,ints);
                    disc = [disc(:); tmp(:)];
                end
                tmp = e.fn(obj.t{n+1},obj.x{n+1},obj.u{n},obj.p,tps,ints);
                disc = [disc(:); tmp(:)];
            end
        end
        
        function rk = rk4q(obj, ocp, q)
            % Integrator parameters
            I0 = casadi.MX.sym('i0'); % Start of integration
            If = casadi.MX.sym('if'); % End of integration
            x0 = casadi.MX.sym('x0', ocp.n_x); % Initial value
            
            % ode parameters
            tt = casadi.MX.sym('t');  % Independent variable
            xx = casadi.MX.sym('x', ocp.n_x); % State variable
            uu = casadi.MX.sym('u', ocp.n_u); % Control input
            pp = casadi.MX.sym('p', ocp.n_p); % Free parameters
            tp = casadi.MX.sym('tp', ocp.n_tp); % Timepoints
            int= casadi.MX.sym('int', ocp.n_int); % integrals
            
            args = {tt,xx,uu,pp,tp,int};
            outs = {ocp.ode(tt,xx,uu,pp), q(tt,xx,uu,pp,tp,int)};
            f = casadi.Function('f', args, outs);
           
            % Expression for RK4 integrator
            tt = I0;
            xx = x0;
            hh = (If-I0)/obj.rk_steps; 
            qq = 0;
            for k=1:obj.rk_steps
                [k1, q1] = f(tt     , xx        , uu, pp, tp, int);
                [k2, q2] = f(tt+hh/2, xx+hh/2*k1, uu, pp, tp, int);
                [k3, q3] = f(tt+hh/2, xx+hh/2*k2, uu, pp, tp, int);
                [k4, q4] = f(tt+hh  , xx+hh  *k3, uu, pp, tp, int);
                tt = tt + hh;
                xx = xx + hh/6*(k1 + 2*k2 + 2*k3 + k4);
                qq = qq + hh/6*(q1 + 2*q2 + 2*q3 + q4);
            end
            rk = casadi.Function('Q', {I0,If,x0,uu,pp,tp,int}, {qq});
        end
        
        function rk = rk4(obj, ocp)            
            % Integrator parameters
            I0 = casadi.MX.sym('I0'); % Start of integration
            If = casadi.MX.sym('If'); % End of integration
            x0 = casadi.MX.sym('x0', ocp.n_x); % Initial value
            
            % ode parameters
            tt  = casadi.MX.sym('t');  % Independent variable
            xx  = casadi.MX.sym('x', ocp.n_x); % State variable
            uu  = casadi.MX.sym('u', ocp.n_u); % Control input
            pp  = casadi.MX.sym('p', ocp.n_p); % Free parameters

            % ode rhs
            f = casadi.Function('f', {tt,xx,uu,pp}, ...
                {ocp.ode(tt,xx,uu,pp)});
           
            % Expression for RK4 integrator
            tt = I0;
            xx = x0;
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
        
        function [t, x, u, p] = vars_at(obj, tp)
                
            if isa(tp, 'yop.ast_independent_initial')
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
        
        function n = n_w(obj)
            n = length(obj.w);
        end
        
    end
end