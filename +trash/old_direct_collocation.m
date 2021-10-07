classdef old_direct_collocation < handle
    properties
        N
        d
        cp
        tau
        
        t0
        tf
        t
        x
        u
        p
        
        J
        g
        h
    end
    
    
    methods
        function obj = old_direct_collocation(ocp, N, d, cp)
            obj.N = N;
            obj.d = d;
            obj.cp = cp;
            obj.tau = [0, casadi.collocation_points(d, cp)];
            
            obj.t0 = casadi.MX.sym('t0');
            obj.tf = casadi.MX.sym('tf');
            obj.t = yop.collocated_time(obj.t0, obj.tf, N);
            obj.x = yop.collocated_mx_state(obj.t0, obj.tf, n_x, N, d, cp);
            obj.u = yop.collocated_mx_control(obj.t0, obj.tf, n_u, N);
            obj.p = casadi.MX.sym('p', n_p);
            
            obj.J = 0;
            obj.g = [];
            obj.h = [];
        end
        
        function sol = solve(obj, bnd)
            nlp = obj.nlp(bnd);
            w0 = ones(size(nlp.w));
            prob = struct('f', nlp.J, 'x', nlp.w, 'g', [nlp.g; nlp.h]);
            solver = casadi.nlpsol('solver', 'ipopt', prob);
            sol = solver( ...
                'x0', w0, ...
                'lbx', nlp.w_lb, ...
                'ubx', nlp.w_ub, ...
                'ubg', [nlp.g_ub; nlp.h_ub], ...
                'lbg', [nlp.g_lb; nlp.h_lb] ...
                ); 
        end
        
        function prob = nlp(obj, bnd)
            [w_lb, w_ub] = obj.set_box_bnd(bnd);
            prob = struct;
            prob.J = obj.J;
            prob.w = obj.w;
            prob.w_ub = w_ub;
            prob.w_lb = w_lb;
            prob.g = obj.g;
            prob.g_ub = zeros(size(obj.g));
            prob.g_lb = zeros(size(obj.g));
            prob.h = obj.h;
            prob.h_ub = zeros(size(obj.h));
            prob.h_lb = -inf(size(obj.h));
        end
        
        function [w_lb, w_ub] = set_box_bnd(obj, bnd)
            t0_lb = bnd.t0_lb;
            t0_ub = bnd.t0_ub;
            tf_lb = bnd.tf_lb;
            tf_ub = bnd.tf_ub;
            
            reps = obj.N*(obj.d + 1) + 1;
            x_lb = repmat(bnd.x_lb, reps, 1);
            x_ub = repmat(bnd.x_ub, reps, 1);
            x_lb(1 : bnd.n_x) = bnd.x0_lb;
            x_ub(1 : bnd.n_x) = bnd.x0_ub;
            x_lb(end - bnd.n_x + 1 : end) = bnd.xf_lb;
            x_ub(end - bnd.n_x + 1 : end) = bnd.xf_ub;
            
            u_lb = repmat(bnd.u_lb, obj.N, 1);
            u_ub = repmat(bnd.u_ub, obj.N, 1);
            u_lb(1 : bnd.n_u) = bnd.u0_lb;
            u_ub(1 : bnd.n_u) = bnd.u0_ub;
            u_lb(end - bnd.n_u + 1 : end) = bnd.uf_lb;
            u_ub(end - bnd.n_u + 1 : end) = bnd.uf_ub;
                        
            p_lb = bnd.p_lb;
            p_ub = bnd.p_ub;
            
            w_ub = vertcat(t0_ub, tf_ub, x_ub, u_ub, p_ub);
            w_lb = vertcat(t0_lb, tf_lb, x_lb, u_lb, p_lb);
        end
        
        function obj = add_ode(obj, ode)
            % Derivative of polynomial
            for n=1:obj.N
                dx = obj.x(n).differentiate();
                for r=2:obj.d+1
                    tt = obj.t(n).evaluate(obj.tau(r));
                    xx = obj.x(n).evaluate(obj.tau(r));
                    uu = obj.u(n).evaluate(obj.tau(r));
                    pp = obj.p;
                    f = ode(tt, xx, uu, pp);
                    dxr = dx.evaluate(obj.tau(r));
                    obj.g = [obj.g; (dxr - obj.dt*f)];
                end
            end
            
            % Continuity
            for n=1:obj.N
                obj.g = [obj.g; ...
                    (obj.x(n).evaluate(1) - obj.x(n+1).evaluate(0))];
            end 
        end
        
        function obj = add_objective( ...
                obj, expr, special_nodes, n_tp, n_int, n_der)
            
            obj.J = obj.parameterize( ...
                expr, ...
                special_nodes, ...
                n_tp, ...
                n_int, ...
                n_der, ...
                obj.tf-obj.t0 ...
                );
        end
        
        function obj = add_eq(obj, expr, special_nodes, n_tp, n_int, n_der)
            disc = obj.parameterize( ...
                expr, ...
                special_nodes, ...
                n_tp, ...
                n_int, ...
                n_der, ...
                obj.tf-obj.t0 ...
                );
            obj.g = [obj.g; disc(:)];
        end
        
        function obj=add_ieq(obj, expr, special_nodes, n_tp, n_int, n_der)
            disc = obj.parameterize( ...
                expr, ...
                special_nodes, ...
                n_tp, ...
                n_int, ...
                n_der, ...
                obj.tf-obj.t0 ...
                );
            obj.h = [obj.h; disc(:)];
        end
        
        function disc = parameterize(obj, expr, sn, n_tp, n_int, n_der, T)
            [tps,ints] = obj.param_special_nodes(sn, n_tp, n_int, n_der,T);
            disc = obj.parameterize_expression(expr, tps, ints);
        end
        
        function [tps,ints]=param_special_nodes(obj,sn,n_tp,n_int,n_der,T)
            tps = [];
            ints = [];
            for node = sn
                tps  = [tps;  zeros(n_tp  - length(tps), 1)];
                ints = [ints; zeros(n_int - length(ints), 1)];
                switch node.type
                    case yop.ocp_expr.tp
                        tp = obj.parameterize_timepoint(node, tps, ints);
                        tps = [tps; tp(:)];
                        
                    case yop.ocp_expr.int
                        int = obj.parameterize_integral(node, tps, ints);
                        ints = [ints; int(:)];
                        
                    case yop.ocp_expr.der
                        error(yop.msg.not_implemented);
                        warning(['When implementing do not forget to ' ...
                            'multiply derivative with step length']);
                    otherwise
                        error(yop.msg.unexpected_error);
                end
            end
        end
        
        function val = parameterize_timepoint(obj, tp, tps, ints, T)
            val = tp.fn( ...
                obj.t.value(tp.timepoint, T), ...
                obj.x.value(tp.timepoint, T), ...
                obj.u.value(tp.timepoint, T), ...
                obj.p, tps, ints);
        end
        
        function I = parameterize_integral(obj, int, tps, ints)
            I = 0;
            for n=1:obj.N
                yval = [];
                for r=1:obj.d+1
                    tt = obj.t(n).evaluate(obj.tau(r));
                    xx = obj.x(n).evaluate(obj.tau(r));
                    uu = obj.u(n).evaluate(obj.tau(r));
                    pp = obj.p;
                    val_r = int.fn(tt, xx, uu, pp, tps, ints);
                    yval = [yval, val_r(:)];
                end
                lp = yop.lagrange_polynomial(obj.tau, yval).integrate();
                I = I + lp.evaluate(1)*obj.dt;
            end
        end
        
        function disc = parameterize_expression(obj, expr, tps, ints)
            if is_transcription_invariant(expr)
                tt = obj.t(1).evaluate(0);
                xx = obj.x(1).evaluate(0);
                uu = obj.u(1).evaluate(0);
                pp = obj.p;
                % dd = der(1).evaluate(0);
                disc = expr.fn(tt, xx, uu, pp, tps, ints);
            else
                disc = [];
                for n=1:obj.N
                    tt = obj.t(n).evaluate(0);
                    xx = obj.x(n).evaluate(0);
                    uu = obj.u(n).evaluate(0);
                    pp = obj.p;
                    % dd = der(n).evaluate(0);
                    disc = [disc, expr.fn(tt, xx, uu, pp, tps, ints)];
                    if expr.is_hard
                        for r=2:length(obj.tau)
                            tt = obj.t(n).evaluate(obj.tau(r));
                            xx = obj.x(n).evaluate(obj.tau(r));
                            uu = obj.u(n).evaluate(obj.tau(r));
                            disc = [disc, ...
                                expr.fn(tt, xx, uu, pp, tps, ints)];
                        end
                    end
                end
                tt = obj.t(obj.N + 1).evaluate(0);
                xx = obj.x(obj.N + 1).evaluate(0);
                uu = obj.u(obj.N + 1).evaluate(0);
                pp = obj.p;
                disc = [disc, expr.fn(tt,xx,uu,pp,tps,ints)];
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
            v = vertcat(obj.t0, obj.tf, vec(obj.x), vec(obj.u), obj.p);
        end
        
        function n = dt(obj)
            n = (obj.tf-obj.t0)/obj.N;
        end
        
        function [t, x, u, p] = ocp_vars(obj, w)
            time = casadi.Function('t', {obj.w}, {mat(obj.t)});
            t = full(time(w));
            
            state = casadi.Function('x', {obj.w}, {mat(obj.x)});
            x = full(state(w));
            
            control = casadi.Function('u', {obj.w}, {mat(obj.u)});
            u = full(control(w));
            
            parameter = casadi.Function('p', {obj.w}, {obj.p});
            p = full(parameter(w))';
        end
        
        function [t0, tf, t, x, u, p] = num_vars(obj, w)
            [tt,xx,uu,p] = obj.ocp_vars(w);
            t0 = tt(1);
            tf = tt(end);
            t = yop.collocated_time(t0, tf, obj.N);
            x = yop.collocated_num_state(t0, tf, xx, obj.N, obj.d, obj.cp);
            u = yop.collocated_num_control(t0, tf, uu, obj.N);
        end
        
    end
    
    
end
