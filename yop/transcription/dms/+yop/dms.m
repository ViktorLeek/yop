classdef dms < handle
    properties
        N
        rk_steps
        ocp
        J
        eq
        ieq
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
        tps
        ints
        diffcon
    end
    methods
        
        function obj = dms(ocp, N, rk_steps)
            obj.ocp = ocp;
            obj.N = N;
            obj.rk_steps = rk_steps;
            obj.build_vars();
            obj.set_box_bnd();
            obj.discretize_ode();
            obj.parameterize_timepoints_and_integrals();
            obj.parameterize_expressions();
        end
        
        function sol = solve(obj)
            sol = obj.solve_nlp();
        end
        
        function sol = solve_nlp(obj)
            w0 = ones(size(obj.w));
            d = obj.diffcon;
            g = obj.ocp.eq.vertcat_disc();
            h = obj.ocp.ieq.vertcat_disc();
            nlp.f = obj.ocp.objective.disc;
            nlp.x = obj.w;
            nlp.g = [d; g; h];
            solver = casadi.nlpsol('solver', 'ipopt', nlp);
            sol = solver( ...
                'x0', w0, ...
                'lbx', obj.w_lb, ...
                'ubx', obj.w_ub, ...
                'ubg', [zeros(size([d;g])); zeros(size(h))], ...
                'lbg', [zeros(size([d;g])); -inf(size(h))] ...
                );
        end
        
        function obj = build_vars(obj)
            % BUILD_VARS - Build nlp variable vector
            
            % Time horizon
            obj.t0 = casadi.MX.sym('t0'); %obj.ocp.t0;
            obj.tf = casadi.MX.sym('tf'); %obj.ocp.tf;
            
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
                obj.x{k} = casadi.MX.sym(['x_' num2str(k)], obj.ocp.n_x);
            end
            
            % Control intervals
            obj.u = cell(obj.N, 1);
            for k=1:obj.N
                obj.u{k} = casadi.MX.sym(['u_' num2str(k)], obj.ocp.n_u);
            end
            
            % Parameters
            obj.p = casadi.MX.sym('p', obj.ocp.n_p);
            
            % Timepoints and integrals
            obj.tps = casadi.MX.sym('tp', obj.ocp.n_tp);
            obj.ints = casadi.MX.sym('int', obj.ocp.n_int);
        end
        
        function obj = set_box_bnd(obj)
            obj.t0_lb = obj.ocp.t0_lb;
            obj.t0_ub = obj.ocp.t0_ub;
            obj.tf_lb = obj.ocp.tf_lb;
            obj.tf_ub = obj.ocp.tf_ub;
            
            obj.x_lb = repmat(obj.ocp.x_lb, obj.N+1, 1);
            obj.x_ub = repmat(obj.ocp.x_ub, obj.N+1, 1);
            obj.x_lb(1:obj.ocp.n_x) = obj.ocp.x0_lb;
            obj.x_ub(1:obj.ocp.n_x) = obj.ocp.x0_ub;
            obj.x_lb(end-obj.ocp.n_x+1:end) = obj.ocp.xf_lb;
            obj.x_ub(end-obj.ocp.n_x+1:end) = obj.ocp.xf_ub;
            
            obj.u_lb = repmat(obj.ocp.u_lb, obj.N, 1);
            obj.u_ub = repmat(obj.ocp.u_ub, obj.N, 1);
            obj.u_lb(1:obj.ocp.n_u) = obj.ocp.u0_lb;
            obj.u_ub(1:obj.ocp.n_u) = obj.ocp.u0_ub;
            obj.u_lb(end-obj.ocp.n_u+1:end) = obj.ocp.uf_lb;
            obj.u_ub(end-obj.ocp.n_u+1:end) = obj.ocp.uf_ub;
                        
            obj.p_lb = obj.ocp.p_lb;
            obj.p_ub = obj.ocp.p_ub;
        end
        
        function obj = discretize_ode(obj)
            F = obj.rk4();
            % Integrate trajectory individually on segments and introduce
            % continuity constraints
            for k=1:(obj.N)
                obj.add_diffcon(obj.x{k+1} - ...
                    F(obj.t{k}, obj.t{k+1}, obj.x{k}, obj.u{k}, obj.p));
            end
        end
        
        function obj = add_diffcon(obj, expr)
            if isempty(obj.diffcon)
                obj.diffcon = expr;
            else
                obj.diffcon = [obj.diffcon(:); expr(:)];
            end
        end
        
        function obj = parameterize_timepoints_and_integrals(obj)
            % Important that this is in topological order, should be
            % handled by ocp. The reason is that they need to be
            % parameterized in order, otherwise there is no gurantee of the
            % results. As more and more timepoints are discretized they can
            % start appearing in expressions with other timepoints, and
            % since this is done in topological order, there is no
            % conflict.
            for tp = obj.ocp.timepoints
                [tt,xx,uu,pp] = obj.vars_at(tp.timepoint);
                tvec = obj.ocp.timepoints.get_disc();
                ivec = obj.ocp.integrals.get_disc();
                tmp = tp.fn(tt,xx,uu,pp,tvec,ivec);
                tp.disc = tmp(:);
            end
            for int = obj.ocp.integrals
                rk = obj.rk4q(int.fn);
                tvec = obj.ocp.timepoints.get_disc();
                ivec = obj.ocp.integrals.get_disc();
                I = 0;
                for n=1:obj.N
                    I = I + rk(obj.t{n}, obj.t{n+1}, obj.x{n}, ...
                        obj.u{n}, obj.p, tvec, ivec);
                end
                int.disc = I(:);
            end
        end
        
        function obj = parameterize_expressions(obj)
            tvec = obj.ocp.timepoints.get_disc();
            ivec = obj.ocp.integrals.get_disc();
            for c = [obj.ocp.objective, obj.ocp.eq, obj.ocp.ieq]
                if all(is_transcription_invariant(c))
                    % Since the expression is transcrition invariant, no of
                    % the variable t,x,u,p will reach the final expression,
                    % so any argument suffices. Instead timepoint,
                    % integrals and constants reach the final expression.
                    tmp = c.fn(obj.t{1},obj.x{1},obj.u{1},obj.p,tvec,ivec);
                    c.disc = tmp(:);
                else
                    grid_val = cell(obj.N+1, 1);
                    for n=1:obj.N
                        tmp = c.fn(obj.t{n}, obj.x{n}, obj.u{n}, obj.p, ... 
                            obj.w, tvec, ivec);
                        grid_val{n} = tmp(:);
                    end
                    tmp = c.fn(obj.t{obj.N+1}, obj.x{obj.N+1}, ...
                        obj.u{obj.N}, obj.p, tvec, ivec);
                    grid_val{obj.N+1} = tmp(:);
                    c.disc = vertcat(grid_val{:});
                end
            end
        end
        
        function rk = rk4q(obj, qfn)            
            nx = obj.ocp.n_x; nu = obj.ocp.n_u; np = obj.ocp.n_p;
            ntp = obj.ocp.n_tp; nint = obj.ocp.n_int;
            
            % Integrator parameters
            I0 = casadi.MX.sym('i0'); % Start of integration
            If = casadi.MX.sym('if'); % End of integration
            x0 = casadi.MX.sym('x0', nx); % Initial value
            
            % ode parameters
            t_ = casadi.MX.sym('t');  % Independent variable
            x_ = casadi.MX.sym('x', nx); % State variable
            u_ = casadi.MX.sym('u', nu); % Control input
            p_ = casadi.MX.sym('p', np); % Free parameters
            tp = casadi.MX.sym('tp', ntp); % Timepoints
            int= casadi.MX.sym('int', nint); % integrals
            
            args = {t_,x_,u_,p_,tp,int};
            outs = {obj.ocp.ode(t_,x_,u_,p_), qfn(t_,x_,u_,p_,tp,int)};
            f = casadi.Function('f', args, outs);
           
            % Expression for RK4 integrator
            tt = I0;
            xx = x0;
            uu = u_;
            pp = p_;
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
        
        function rk = rk4(obj)            
            nx = obj.ocp.n_x; nu = obj.ocp.n_u; np = obj.ocp.n_p;
            
            % Integrator parameters
            I0 = casadi.MX.sym('I0'); % Start of integration
            If = casadi.MX.sym('If'); % End of integration
            x0 = casadi.MX.sym('x0', nx); % Initial value
            
            % ode parameters
            % (Awkward naming, but otherwise Matlab's linter kills my eyes)
            t_  = casadi.MX.sym('t');  % Independent variable
            x_  = casadi.MX.sym('x', nx); % State variable
            u_  = casadi.MX.sym('u', nu); % Control input
            p_  = casadi.MX.sym('p', np); % Free parameters

            % ode rhs
            f = casadi.Function('f', {t_,x_,u_,p_}, ...
                {obj.ocp.ode(t_,x_,u_,p_)});
           
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
        
        function bd = w_lb(obj)
            bd = vertcat(obj.t0_lb, obj.tf_lb, obj.x_lb, obj.u_lb, obj.p_lb);
        end
        
        function bd = w_ub(obj)
            bd = vertcat(obj.t0_ub, obj.tf_ub, obj.x_ub, obj.u_ub, obj.p_ub);
        end
        
        function n = n_w(obj)
            n = length(obj.w);
        end
        
    end
end