classdef direct_multiple_shooting < handle
    properties
        ivals
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
            obj.ivals = N;
            obj.rk_steps = rk_steps;
        end
        
        function [t, x, u, p] = to_numeric(obj, w)
            time = casadi.Function('x', {obj.w}, {vertcat(obj.t{:})});
            t = full(time(w));
            
            state = casadi.Function('x', {obj.w}, {horzcat(obj.x{:})});
            x = full(state(w))';
            
            control = casadi.Function('x', {obj.w}, {horzcat(obj.u{:})});
            u = full(control(w))';
            
            parameter = casadi.Function('x', {obj.w}, {horzcat(obj.p{:})});
            p = full(parameter(w))';
        end
        
        function nlp = transcribe(obj, ocp)
            obj.build_vars(ocp);
            [w_lb, w_ub] = obj.set_box_bnd(ocp);
            
            cc = obj.discretize_ode(ocp); % Continuity constraints
            [tps, ints] = parameterize_special_nodes(obj, ocp);
            J = obj.parameterize_objective(ocp, tps, ints);
            g = obj.parameterize_equality(ocp, tps, ints);
            h = obj.parameterize_inequality(ocp, tps, ints);
            
            nlp = struct;
            nlp.J = J;
            nlp.w = obj.w;
            nlp.w_ub = w_ub;
            nlp.w_lb = w_lb;
            nlp.g = [cc; g];
            nlp.g_ub = zeros(size([cc; g]));
            nlp.g_lb = zeros(size([cc; g]));
            nlp.h = h;
            nlp.h_ub = zeros(size(h));
            nlp.h_lb = -inf(size(h));
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
                cc = [cc; ck(:)];
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
            rk = obj.rk4q(ocp, int.fn);
            I = 0;
            for n=1:obj.N
                ti  = obj.t{n};
                tii = obj.t{n+1};
                x0  = obj.x{n};
                uu  = obj.u{n};
                pp  = obj.p;
                I = I + rk(ti, tii, x0, uu, pp, tp_tmp, int_tmp);
            end
            ints = [ints; I(:)];
            
        end
        
%         function tps = parameterize_timepoints(obj, ocp)
%             % Important that this is in topological order, should be
%             % handled by ocp. The reason is that they need to be
%             % parameterized in order, otherwise there is no gurantee of the
%             % results. As more and more timepoints are discretized they can
%             % start appearing in expressions with other timepoints, and
%             % since this is done in topological order, there is no
%             % conflict.
%             ints = zeros(ocp.n_int, 1);
%             tps = [];
%             for k = 1:length(ocp.timepoints)
%                 tp_k = ocp.timepoints(k);
%                 [tt, xx, uu, pp] = obj.vars_at(tp_k.timepoint, ocp);
%                 tmp = [tps; zeros(ocp.n_tp - length(tps), 1)];
%                 tp_k = tp_k.fn(tt, xx, uu, pp, tmp, ints);
%                 tps = [tps; tp_k(:)];
%             end
%         end
%         
%         function ints = parameterize_integrals(obj, ocp, tps)
%             % Important that this is in topological order, should be
%             % handled by ocp. The reason is that they need to be
%             % parameterized in order, otherwise there is no gurantee of the
%             % results. As more and more timepoints are discretized they can
%             % start appearing in expressions with other timepoints, and
%             % since this is done in topological order, there is no
%             % conflict.
%             ints = [];
%             for k = 1:length(ocp.integrals)
%                 rk = obj.rk4q(ocp, ocp.integrals(k).fn);
%                 tmp = [ints; zeros(ocp.n_int - length(ints), 1)];
%                 I = 0;
%                 for n=1:obj.N
%                     ti  = obj.t{n};
%                     tii = obj.t{n+1};
%                     x0  = obj.x{n};
%                     uu  = obj.u{n};
%                     pp  = obj.p;
%                     I = I + rk(ti, tii, x0, uu, pp, tps, tmp);
%                 end
%                 ints = [ints; I(:)];
%             end
%         end
        
        function J = parameterize_objective(obj, ocp, tps, ints)
            J = obj.parameterize_expression(ocp.objective, tps, ints);
        end
        
        function g = parameterize_equality(obj, ocp, tps, ints)
            g = [];
            for k=1:length(ocp.equality_constraints)
                g_k = obj.parameterize_expression( ...
                    ocp.equality_constraints(k), tps, ints);
                g = [g; g_k(:)];
            end
        end
        
        function h = parameterize_inequality(obj, ocp, tps, ints)
            h = [];
            for k=1:length(ocp.inequality_constraints)
                h_k = obj.parameterize_expression( ...
                    ocp.inequality_constraints(k), tps, ints);
                h = [h; h_k(:)];
            end
        end
        
        function disc = parameterize_expression(obj, expr, tps, ints)
            if is_transcription_invariant(expr)
                tt = obj.t{1};
                xx = obj.x{1};
                uu = obj.u{1};
                pp = obj.p;
                disc = expr.fn(tt, xx, uu, pp, tps, ints);
            else
                disc = [];
                for n=1:obj.N
                    tt = obj.t{n};
                    xx = obj.x{n};
                    uu = obj.u{n};
                    pp = obj.p;
                    disc = [disc; expr.fn(tt, xx, uu, pp, tps, ints)];
                end
                tt = obj.t{obj.N+1};
                xx = obj.x{obj.N+1};
                uu = obj.u{obj.N};
                pp = obj.p;
                disc = [disc; expr.fn(tt, xx, uu, pp, tps, ints)];
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
        
        function [t, x, u, p] = vars_at(obj, tp, ocp)
                
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
            
            [fixed, T] = fixed_horizon(ocp);            
            
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
        
        function n = N(obj)
            n = obj.ivals;
        end
    end
end