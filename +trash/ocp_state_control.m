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
        
        function v = vec(obj)
            v = [];
            for o=obj
                v = [v; o.ast];
            end
        end
        
        function w = weight(obj)
            w = [];
            for o=obj
                w = [w; o.ast.weight];
            end
        end
        
        function os = offset(obj)
            os = [];
            for o=obj
                os = [os; o.ast.offset];
            end
        end
        
        function ID = ids(obj)
            ID = [];
            for o = obj
                ID = [ID, o.ast.id];
            end
        end
    end
end