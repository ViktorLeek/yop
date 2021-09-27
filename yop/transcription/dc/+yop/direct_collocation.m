classdef direct_collocation < handle
    properties
        ivals
        points
        polydeg
        tau
        t0
        tf
        t
        x
        u
        p
    end
    methods
        function obj = direct_collocation(ivals, polydeg, points, n_x, n_u, n_p)
            obj.ivals = ivals;
            obj.points = points;
            obj.polydeg = polydeg;
            obj.build_vars(n_x, n_u, n_p);
        end
        
        function [t, x, u, p, tx] = to_numeric(obj, w)
            tt = [];
            for n=1:obj.N
                for r=1:obj.d+1
                    tt = [tt(:); obj.t{n,r}];
                end
            end
            tt = [tt(:); obj.t{obj.N+1,1}];
            time = casadi.Function('x', {obj.w}, {tt});
            tx = full(time(w));
            
            tt = [];
            for n=1:obj.N+1
                tt = [tt; obj.t{n,1}];
            end
            time = casadi.Function('x', {obj.w}, {vertcat(tt)});
            t = full(time(w));
            
            xx = [];
            for xk=obj.x
                xx = [xx; xk.y'];
            end
            state = casadi.Function('x', {obj.w}, {xx});
            x = full(state(w));
            
            control = casadi.Function('x', {obj.w}, {horzcat(obj.u{:})});
            u = full(control(w))';
            
            parameter = casadi.Function('x', {obj.w}, {horzcat(obj.p{:})});
            p = full(parameter(w))';
        end
        
        function nlp = transcribe(obj, ocp)
            [w_lb, w_ub] = obj.set_box_bnd(ocp);
            
            c = obj.discretize_ode(ocp);
            [tps, ints] = parameterize_special_nodes(obj, ocp);
            J = obj.parameterize_objective(ocp, tps, ints);
            g = obj.parameterize_equality(ocp, tps, ints);
            h = obj.parameterize_inequality(ocp, tps, ints);
            
            nlp = struct;
            nlp.J = J;
            nlp.w = obj.w;
            nlp.w_ub = w_ub;
            nlp.w_lb = w_lb;
            nlp.g = [c; g];
            nlp.g_ub = zeros(size([c; g]));
            nlp.g_lb = zeros(size([c; g]));
            nlp.h = h;
            nlp.h_ub = zeros(size(h));
            nlp.h_lb = -inf(size(h));
        end
        
        function obj = build_vars(obj, n_x, n_u, n_p)
           obj.tau = [0, ...
               casadi.collocation_points(obj.polydeg, obj.points)];
           
           % Independent
           obj.t0 = casadi.MX.sym('t0');
           obj.tf = casadi.MX.sym('tf');
           obj.t = cell(obj.N+1, obj.d+1);
           for n=1:obj.N
               for r=1:(obj.d+1)
                   obj.t{n,r} = obj.t0 + obj.h*(n-1) + obj.h*obj.tau(r);
               end
           end
           obj.t{obj.N+1,1} = obj.tf;
           
           % State
           obj.x = yop.collocated_state( ...
               'x', n_x, obj.N, obj.d, obj.points);
           
           % Control
           obj.u = cell(obj.N, 1);
           for n=1:obj.N
               obj.u{n} = casadi.MX.sym(['u_', num2str(n)], n_u);
           end
           
           % Parameter
           obj.p = casadi.MX.sym('p', n_p);
        end
        
        function [w_lb, w_ub] = set_box_bnd(obj, ocp)
            t0_lb = ocp.t0_lb;
            t0_ub = ocp.t0_ub;
            tf_lb = ocp.tf_lb;
            tf_ub = ocp.tf_ub;
            
            reps = obj.N*(obj.d + 1) + 1;
            x_lb = repmat(ocp.x_lb, reps, 1);
            x_ub = repmat(ocp.x_ub, reps, 1);
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
        
        function c = discretize_ode(obj, ocp)
            % Derivative of polynomial
            c = [];
            for n=1:obj.N
                dx = obj.x(n).differentiate();
                for r=2:obj.d+1
                    tt = obj.t{n,r};
                    xx = obj.x(n).evaluate(obj.tau(r));
                    uu = obj.u{n};
                    pp = obj.p;
                    f = ocp.ode(tt, xx, uu, pp);
                    dxr = dx.evaluate(obj.tau(r));
                    c = [c; (dxr - obj.h*f)];
                end
            end
            
            % Continuity
            for n=1:obj.N
                c = [c; (obj.x(n).evaluate(1) - obj.x(n+1).evaluate(0))];
            end 
            
        end
        
        function [tps, ints] = parameterize_special_nodes(obj, ocp)
            tps = [];
            ints = [];
            for node = ocp.special_nodes
                switch node.type
                    case yop.ocp_expr.tp
                        tps = obj.parameterize_timepoint( ...
                            node, tps, ints, ocp);
                    case yop.ocp_expr.int
                        ints = obj.parameterize_integral( ...
                            node, tps, ints, ocp);
                    case yop.ocp_expr.der
                        error(yop.msg.not_implemented);
                        warning(['When implementing dont forget to ' ...
                            'multiply derivative with step length']);
                    otherwise
                        error(yop.msg.unexpected_error);
                end
            end
        end
        
        function tps = parameterize_timepoint(obj, tp, tps, ints, ocp)
            tp_tmp  = [tps;  zeros(ocp.n_tp  - length(tps), 1)];
            int_tmp = [ints; zeros(ocp.n_int - length(ints), 1)];
            [tt, xx, uu, pp] = obj.vars_at(tp.timepoint, ocp);
            value = tp.fn(tt, xx, uu, pp, tp_tmp, int_tmp);
            tps = [tps; value(:)];
        end
        
        function ints = parameterize_integral(obj, int, tps, ints, ocp)
            tp_tmp  = [tps;  zeros(ocp.n_tp  - length(tps), 1)];
            int_tmp = [ints; zeros(ocp.n_int - length(ints), 1)];
            I = 0;
            for n=1:obj.N
                yval = [];
                for r=1:obj.d+1
                    tt = obj.t{n, r};
                    xx = obj.x(n).evaluate(obj.tau(r));
                    uu = obj.u{n};
                    pp = obj.p;
                    val_r = int.fn(tt, xx, uu, pp, tp_tmp, int_tmp);
                    yval = [yval, val_r(:)];
                end
                lp = yop.lagrange_polynomial(obj.tau, yval).integrate();
                % Important detail is to multiply with step length
                I = I + lp.evaluate(1)*obj.h;
            end
            ints = [ints; I(:)];
        end
        
        function J = parameterize_objective(obj, ocp, tps, ints)
            J = obj.parameterize_expression(ocp.objective, tps, ints);
        end
        
        function g = parameterize_equality(obj, ocp, tps, ints)
            g = [];
            for eq_k=1:length(ocp.equality_constraints)
                g_k = obj.parameterize_expression(eq_k, tps, ints);
                g = [g; g_k(:)];
            end
        end
        
        function h = parameterize_inequality(obj, ocp, tps, ints)
            h = [];
            for ieq_k=ocp.inequality_constraints
                h_k = obj.parameterize_expression(ieq_k, tps, ints);
                h = [h; h_k(:)];
            end
        end
        
        function disc = parameterize_expression(obj, expr, tps, ints)
            if is_transcription_invariant(expr)
                tt = obj.t{1, 1};
                xx = obj.x(1).evaluate(0);
                uu = obj.u{1};
                pp = obj.p;
                disc = expr.fn(tt, xx, uu, pp, tps, ints);
            else
                % hard
                disc = [];
                for n=1:obj.N
                    tt = obj.t{n,1};
                    xx = obj.x(n).evaluate(0);
                    uu = obj.u{n};
                    pp = obj.p;
                    disc = [disc; expr.fn(tt, xx, uu, pp, tvec, ivec)];
                    if is_hard(expr)
                        for r=2:length(obj.tau)
                            tt = obj.t{n,r};
                            xx = obj.x(n).evaluate(obj.tau(r));
                            disc = [disc; ...
                                expr.fn(tt, xx, uu, pp, tvec, ivec)];
                        end
                    end
                end
                tt = obj.t{obj.N + 1, 1};
                xx = obj.x(obj.N + 1).evaluate(0);
                uu = obj.u{obj.N};
                pp = obj.p;
                disc = [disc; c.fn(tt,xx,uu,pp,tvec,ivec)];
            end
        end
        
        function [t, x, u, p] = vars_at(obj, tp, ocp)
                
            if isa(tp, 'yop.ast_independent_initial')
                t = obj.t0;
                x = obj.x(1).y(:,1); % Need an abstraction here
                u = obj.u{1};
                p = obj.p;
                return;
            elseif isa(tp, 'yop.ast_independent_final')
                t = obj.tf;
                x = obj.x(end).y; % y-value at N+1 is a vector 
                u = obj.u{end};
                p = obj.p;
                return;
            end
            
            [fixed, T] = fixed_horizon(ocp);            
            
            if ~fixed
                error('[Yop] Problem does not have a fixed horizon')
            end
            
            dt = T/obj.N; % grid step_length    
            if yop.EQ(rem(tp,dt), 0, 1e-3) % 1e-3 should be a setting (also for DMS)
                % grid point
                n = round(tp/dt);
                t = obj.t{n};
                x = obj.x(n).y(:,1);
                u = obj.u{min(obj.N, n)};
                p = obj.p;
                
            else
                % With the first statement and this floor operation n!=N+1,
                % unless t is outside the horizon, in which case an error
                % is desirable.
                % Integrate up to timepoint
                n = floor(tp/dt);
                x = obj.x(n).evaluate(tp);
                u = obj.u{n};
                p = obj.p;
            end
        end
        
        function obj = add_diffcon(obj, expr)
            if isempty(obj.diffcon)
                obj.diffcon = expr(:);
            else
                obj.diffcon = [obj.diffcon; expr(:)];
            end
        end
        
        function v = w(obj)
            xx = [];
            for xk=obj.x
                xx = [xx; xk.y(:)];
            end
            v = vertcat(obj.t0, obj.tf, xx, obj.u{:}, obj.p);
        end
        
        function n = N(obj)
            n = obj.ivals;
        end
        
        function n = d(obj)
            n = obj.polydeg;
        end
        
        function n = h(obj)
            n = (obj.tf-obj.t0)/obj.N;
        end
        
    end
end
