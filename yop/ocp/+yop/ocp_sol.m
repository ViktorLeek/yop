classdef ocp_sol < handle
    properties
        sol
        tro
        
        t
        x
        z
        u
        p
    end
    methods
        function obj = ocp_sol(sol, tro, t, x, z, u, p)
            obj.sol = sol;
            obj.tro = tro;
            obj.t = t;
            obj.x = x;
            obj.z = z;
            obj.u = u;
            obj.p = p;
        end
        
        function obj = reinit(obj,t,x,z,u,p)
        end
        
        function v = value(obj, expr)
            e = yop.ocp_expr(expr);
            
            [sn, n_tp, n_int, n_der] = ...
                mx_functionx(e, obj.t, obj.x, obj.u, obj.p);
            
            [tps, ints] = ...
                obj.tro.parameterize_special_nodes(sn, n_tp, n_int, n_der);
            
            tpfn = casadi.Function('v', {obj.tro.w}, {tps});
            tpv = full(tpfn(obj.sol.x));
            
            intfn = casadi.Function('v', {obj.tro.w}, {ints});
            intv = full(intfn(obj.sol.x));
            
            disc = obj.tro.parameterize_expression(e, tpv, intv);
            fn = casadi.Function('F', {obj.tro.w}, {disc});
            v = full(fn(obj.sol.x));
            
        end
        
        function varargout = plot(obj, varargin)
        end
    end
end