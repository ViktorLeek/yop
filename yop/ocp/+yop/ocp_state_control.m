classdef ocp_state_control < handle
    properties
        ast
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
        function obj = ocp_state_control(ast)
            obj.ast = ast;
            obj.mx = yop.cx(['ocp_', ast.name]);
            obj.sym = sym(ast.name, size(ast));
        end
        
        function obj = set_value(obj, value)
            obj.ast.m_value = value;
        end
        
        
        function obj = set_sym(obj)
            for o=obj
                o.ast.m_value = o.sym;
            end
        end
        
        function obj = set_mx(obj)
            for o=obj
                o.ast.m_value = o.mx;
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
                vec = [vec; o.ast];
            end
        end
    end
end