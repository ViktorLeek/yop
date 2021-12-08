classdef ocp_var < handle
    properties
        var
        mx
        sym
        ub0  % initial upper bound
        lb0  % initial lower bound
        ub   % upper bound
        lb   % lower bound
        ubf  % final upper bound
        lbf  % final lower bound
    end
    methods
        function obj = ocp_var(ast)       
            obj.var = ast;
            obj.mx = yop.cx(['ocp_', ast.name]);
            obj.sym = sym(ast.name, size(ast));
        end
        
        function n = n_elem(obj)
            n = length(obj);
        end
        
        function obj = set_value(obj, value)
            obj.var.m_value = value;
        end
        
        
        function obj = set_sym(obj)
            for o=obj
                o.var.m_value = o.sym;
            end
        end
        
        function obj = set_mx(obj)
            for o=obj
                o.var.m_value = o.mx;
            end
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for o=obj
                vec = [vec; o.mx];
            end
        end
        
        function vec = vec(obj)
            vec = [];
            for o=obj
                vec = [vec; o.var];
            end
        end
        
    end
end