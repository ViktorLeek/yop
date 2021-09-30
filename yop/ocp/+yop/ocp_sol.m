classdef ocp_sol < handle
    properties
        sol
        solver
        nlp
        
        ocp_t
        ocp_x
        ocp_z
        ocp_u
        ocp_p
    end
    methods
        function obj = ocp_sol(sol, solver, nlp, t, x, z, u, p)
            obj.sol = sol;
            obj.solver = solver;
            obj.nlp = nlp;
            obj.ocp_t = t;
            obj.ocp_x = x;
            obj.ocp_z = z;
            obj.ocp_u = u;
            obj.ocp_p = p;
        end
        
        function obj = reinit(obj,t,x,z,u,p)
        end
        
        function v = value(obj, expr)
            % This discretizes the expression in the same way that the
            % transcription method does it. The idea is to see what the
            % solver sees.
            
            expr = yop.ocp_expr(expr);
            
            [vars,tps,ints,ders,sn] = yop.ocp.find_special_nodes(expr.ast);
            
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    vars{k}.m_value = obj.ocp_t.mx;
                end
            end
            
            args = { ...
                mx_vec(obj.ocp_t), ...
                mx_vec(obj.ocp_x), ...
                mx_vec(obj.ocp_u), ...
                mx_vec(obj.ocp_p), ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx([obj.ocp_t, obj.ocp_x, obj.ocp_u, obj.ocp_p]);
            set_mx([tps, ints, ders]);
            
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            
            expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});
            
            t0 = full(obj.sol.x(1));
            tf = full(obj.sol.x(2));
            [tpv, intv] = yop.param_special_nodes(sn, n_elem(tps), ...
                n_elem(ints), obj.nlp.N, obj.nlp.tau, obj.nlp.dt, t0, ...
                tf, obj.nlp.t, obj.nlp.x, obj.nlp.u, obj.nlp.p);
            
            disc = yop.parameterize_expression(expr, obj.nlp.N, ...
                obj.nlp.tau, obj.nlp.t, obj.nlp.x, obj.nlp.u, ...
                obj.nlp.p, tpv, intv, []);
            
            f = casadi.Function('f', {obj.nlp.w}, {disc});
            v = full(f(obj.sol.x));
        end
        
        function v = valuex(obj, expr, refinement)
            expr = yop.ocp_expr(expr);
            
            %% This could be part of the contructor and tt, xx ... could be
            %  instance variables
            time = casadi.Function('t', {obj.nlp.w}, {mat(obj.nlp.t)});
            tt = full(time(obj.sol.x));
            
            state = casadi.Function('x', {obj.nlp.w}, {mat(obj.nlp.x)});
            xx = full(state(obj.sol.x));
            
            control = casadi.Function('u', {obj.nlp.w}, {mat(obj.nlp.u)});
            uu = full(control(obj.sol.x));
            
            parameter = casadi.Function('p', {obj.nlp.w}, {obj.nlp.p});
            pp = full(parameter(obj.sol.x));
            
            tlp = yop.collocated_time(tt(1), tt(end), obj.nlp.N);
            
            y = cell(obj.nlp.N+1,1);
            for n=1:obj.nlp.N
                y{n} = xx(:, (obj.nlp.d+1)*n-obj.nlp.d:(obj.nlp.d+1)*n);
            end
            y{obj.nlp.N+1} = xx(:,end);
            xlp = yop.collocated_expression(obj.nlp.N, obj.nlp.tau, y);
            
            y = cell(obj.nlp.N+1,1);
            for n=1:obj.nlp.N
                y{n} = uu(:, n);
            end
            y{obj.nlp.N+1} = uu(:,end);
            ulp = yop.collocated_expression(obj.nlp.N, 0, y);
            %%
            
            [vars,tps,ints,ders,sn] = yop.ocp.find_special_nodes(expr.ast);
            for k=1:length(vars)
                if isa(vars{k}, 'yop.ast_independent')
                    vars{k}.m_value = obj.ocp_t.mx;
                end
            end
            
            args = { ...
                mx_vec(obj.ocp_t), ...
                mx_vec(obj.ocp_x), ...
                mx_vec(obj.ocp_u), ...
                mx_vec(obj.ocp_p), ...
                mx_vec(tps), ...
                mx_vec(ints), ...
                mx_vec(ders) ...
                };
            
            set_mx([obj.ocp_t, obj.ocp_x, obj.ocp_u, obj.ocp_p]);
            set_mx([tps, ints, ders]);
            
            for node = [tps, ints, ders]
                mx_expr = fw_eval(node.ast.expr);
                node.fn = casadi.Function('fn', args, {mx_expr});
            end
            
            expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});
            
            t0 = full(obj.sol.x(1));
            tf = full(obj.sol.x(2));
            [tpv, intv] = yop.param_special_nodes(sn, n_elem(tps), ...
                n_elem(ints), obj.nlp.N, obj.nlp.tau, obj.nlp.dt, t0, ...
                tf, obj.nlp.t, obj.nlp.x, obj.nlp.u, obj.nlp.p);
            
            tpf = casadi.Function('f', {obj.nlp.w}, {tpv});
            intf = casadi.Function('f', {obj.nlp.w}, {intv});
            
            TP = tpf(obj.sol.x);
            I = intf(obj.sol.x);
            
            if expr.is_transcription_invariant
                v = expr.fn( ...
                    tlp(1).evaluate(0), ...
                    xlp(1).evaluate(0), ...
                    ulp(1).evaluate(0), ...
                    pp, TP, I, []);
            else
                t = tt(1:2:end);
                dt = t(2)-t(1);
                tvec = t;
                for k=1:refinement-1
                    tvec = [tvec; t+k*dt/refinement];
                end
                v = [];
                for tp=tvec(:)'
                    vk = expr.fn( ...
                        tlp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
                        xlp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
                        ulp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
                        pp, TP, I, []);
                    v = [v, vk];
                end
            end
            v = full(v);
        end
        
        function varargout = plot(obj, varargin)
            
            args = {};
            refine = false;
            refinement = 0;
            k = 1;
            while k <= length(varargin)
                if strcmp(varargin{k}, 'refine')
                    refine = true;
                    refinement = varargin{k+1};
                    k = k + 2;
                else
                    args{end+1} = varargin{k};
                    k = k + 1;
                end
            end
            
            for k=1:length(args)
                if isa(args{k}, 'yop.node')
                    if refine
                        args{k} = obj.valuex(varargin{k}, refinement);
                    else
                        args{k} = obj.value(varargin{k});
                    end
                end
            end
            
            h = plot(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
        
        function varargout = stairs(obj, varargin)
            args = {};
            refine = false;
            refinement = 0;
            k = 1;
            while k <= length(varargin)
                if strcmp(varargin{k}, 'refine')
                    refine = true;
                    refinement = varargin{k+1};
                    k = k + 2;
                else
                    args{end+1} = varargin{k};
                    k = k + 1;
                end
            end
            
            for k=1:length(args)
                if isa(args{k}, 'yop.node')
                    if refine
                        args{k} = obj.valuex(varargin{k}, refinement);
                    else
                        args{k} = obj.value(varargin{k});
                    end
                end
            end
            
            h = stairs(args{:});
            if nargout > 0
                varargout{1} = h;
            end 
        end
    end
end