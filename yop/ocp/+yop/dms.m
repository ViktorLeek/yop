classdef dms < handle
    properties
        N
        rk_steps
        ocp 
        ocp_var
        nx
        nu
        np
        nlp
        nlp_var
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
    end
    methods
        function obj = dms(ocp, N, rk_steps)
            obj.ocp = ocp;
            obj.N = N;
            obj.rk_steps = rk_steps;
        end
        
        function [t, x, u, p] = vars_at(obj, tp)
            
            [fixed, T] = fixed_horizon(obj);
            
            if ~fixed
                t=[]; x=[]; u=[]; p=[];
                return;
            end
            
            h = T/obj.N; % grid step_length
            
            if isa(tp, 'yop.independent_initial')
                t = obj.t{1};
                x = obj.x{1};
                u = obj.u{1};
                p = obj.p;
                
            elseif isa(tp, 'yop.ast_independent_final')
                t = obj.t{end};
                x = obj.x{end};
                u = obj.u{end};
                p = obj.p;
            
            elseif abs(rem(tp, h)) <= 1e-9
                % No integration required as the point is close enough to a
                % grid point.
                n = round(tp/h);
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
                n = floor(tp/h); 
                x = ode_int(obj.t{n}, tp, obj.x{n}, obj.u{n}, obj.p{n});
                u = obj.u{n};
                p = obj.p;
            end
        end
        
        function obj = set_nlp_vars(obj)
            % Time horizon
            obj.t0 = obj.ocp.t0;
            obj.tf = obj.ocp.tf;
            
            % grid points
            h = (obj.tf - obj.t0)/obj.N;
            obj.t = cell(obj.N+1, 1);
            obj.t{1} = obj.t0;
            obj.t{end} = obj.tf;
            for k=2:obj.N
                obj.t{k} = obj.t{k-1} + h;
            end
            
            % State at grid points
            obj.x = cell(obj.N+1,1);
            for k=1:(obj.N+1)
                obj.x{k} = casadi.MX.sym(['x_' num2str(k)], obj.nx);
            end
            
            % Control intervals
            obj.u = cell(obj.N, 1);
            for k=1:obj.N
                obj.u{k} = casadi.MX.sym(['u_' num2str(k)], obj.nu);
            end
            
            % Parameters
            obj.p = obj.ocp.p;
        end
        
        function obj = set_box(obj)
            [t0lb, t0ub, tflb, tfub] = obj.ocp.independent_bound();
            obj.t0_lb = t0lb;
            obj.t0_ub = t0ub;
            obj.tf_lb = tflb;
            obj.tf_ub = tfub;
            
            [x0lb, x0ub, xlb, xub, xflb, xfub] = obj.ocp.state_bound();
            obj.x_lb = repmat(xlb, obj.N+1, 1);
            obj.x_ub = repmat(xub, obj.N+1, 1);
            obj.x_lb(1:obj.nx) = x0lb;
            obj.x_ub(1:obj.nx) = x0ub;
            obj.x_lb(end-obj.nx+1:end) = xflb;
            obj.x_ub(end-obj.nx+1:end) = xfub;
            
            [u0lb, u0ub, ulb, uub, uflb, ufub] = obj.ocp.control_bound();
            obj.u_lb = repmat(ulb, obj.N, 1);
            obj.u_ub = repmat(uub, obj.N, 1);
            obj.u_lb(1:obj.nu) = u0lb;
            obj.u_ub(1:obj.nu) = u0ub;
            obj.u_lb(end-obj.nu+1:end) = uflb;
            obj.u_ub(end-obj.nu+1:end) = ufub;
            
            [plb, pub] = obj.ocp.parameter_bound();
            obj.p_lb = plb;
            obj.p_ub = pub; 
        end
        
        function v = w(obj)
            v = vertcat(obj.t0, obj.tf, obj.x{:}, obj.u{:}, obj.p);
        end
        
        function val = integrate_expr(obj, fn)
            int = yop.rk4_dms_expr(fn, obj.ocp.nx, obj.ocp.nu, ...
                obj.ocp.np, obj.length(obj.w), obj.rk_steps);
            val = 0;
            for k=1:obj.N
                val = val + int(...
                    obj.t{k}, obj.t{k+1}, obj.x{k}, obj.u{k}, obj.p,obj.w);
            end
        end
        
        function param_expr(obj, expr)
            [tsort, n_elem] = topological_sort(expr);
            
            % Reset all predecessors so that only nodes in expr is included
            for k=1:n_elem
                reset_pred(tsort{k});
            end
            
            % Set predecessors
            for k=1:n_elem
                set_pred(tsort{k});
            end
            
            % Compute dominators
            for k=n_elem:-1:1
                comp_dom(tsort{k});
            end
            
        end
        
        function e = parameterize_expression(obj, expr)
            % PARAMETERIZE_EXPRESSION
            % Transcribe an ocp expression from the continuous time domain
            % of the ocp to the static domain of the NLP.
            %
            % The implementation is based on a toplogical sort of all the
            % nodes in the expression. Because the sort is topological it
            % only depends on nodes that has been processed earlier or
            % itself, which makes it possible to transcribe the nodes as
            % they are processed
            
            % 1. First make a topological sort
            [topsort, n_elem] = topological_sort(expr);
            
            % 2. Parameterize static parts using the nlp variables at that
            % time instant, integrals using the nlp variables for the
            % entire horizon, and simple time dependent variables using 
            args = {obj.ocp.t, obj.ocp.x, obj.ocp.u, obj.ocp.p, obj.w};
            
            obj.ocp.set_sym();
            for k=1:n_elem
                if isa(topsort{k}, 'yop.ast_timepoint')
                    % Parametrization of timepoint. 
                    % To parameterize it with the values at the timepoint,
                    % the ocp first generates an expression function for
                    % the timepoint expression. The function is evaluated
                    % with values at the timepoint, and the node value is
                    % set to that value. If we don't call forward on this
                    % node (m_value is not set) 
                    fn = obj.ocp.expr_fn(topsort{k}.expr);
                    [tt, xx, uu, pp] = obj.vars_at(topsort{k}.timepoint);
                    topsort{k}.m_value = fn(tt, xx, uu, pp);
                    
                elseif isa(topsort{k}, 'yop.ast_int')
                    integrand = forward(obj);
                    
                    fn = casadi.Function('I', args, {integrand});
                    topsort{k}.m_value = obj.integrate_expr(fn);
                    
                else
                    forward(topsport{k});
                    
                end
            end
            obj.ocp.reset_variables();
            
            
            
            if istimeindependent
                e = {topsort{k}.m_value};
            else
                e = cell(obj.N+1, 1);
                fn = casadi.Function('e', args, {e});
                for k=1:(obj.N)
                    e{k} = ...
                        fn(obj.t{k}, obj.x{k}, obj.u{k}, obj.p{k}, obj.w);
                end
                e{obj.N+1} = fn(...
                    obj.t{obj.N+1}, ...
                    obj.x{obj.N+1}, ...
                    obj.u{obj.N}, ...
                    obj.p{obj.N+1}, ...
                    obj.w);
            end
            
        end
    end
end