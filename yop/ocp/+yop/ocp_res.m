classdef ocp_res < handle
    properties
        independent
        independent_initial
        independent_final
        states
        algebraics
        controls
        parameters
        
        transcriptor
    end
    methods
        function obj = ocp_res()
        end
        
        function obj = reinit(obj,t,x,z,u,p)
        end
        
        function v = t(obj)
        end
        
        function v = x(obj)
        end
        
        function v = u(obj)
        end
        
        function v = z(obj)
        end
        
        function v = p(obj)
        end
        
        function v = value(obj, expr)
            [sorted, tps, ints, ders] = yop.find_special_nodes(expr);
            yop.set_mx_functions(sorted, t, x, z, u, p, tps, ints, ders);
            
            set_mx(obj.variables);
            set_mx(tps);
            set_mx(ints);
            set_mx(ders);
            mx_expr = fw_eval(expr);
            
            fn = casadi.Function('expr', ...
                {t,x,z,u,p,tps,ints,ders},{mx_expr});
            v = obj.transcriptor.value(fn, tps, ints, ders);
        end
        
        function varargout = plot(obj, varargin)
        end
    end
end