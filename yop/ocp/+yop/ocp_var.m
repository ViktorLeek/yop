classdef ocp_var < handle
    properties
        ast
        ub0  % initial upper bound
        lb0  % initial lower bound
        ub   % upper bound
        lb   % lower bound
        ubf  % final upper bound
        lbf  % final lower bound
    end
    
    methods
        function obj = ocp_var(ast)
            obj.ast = ast;
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for o=obj
                vec = [vec; o.ast.m_value];
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