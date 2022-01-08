classdef ivp_var < handle
    properties
        ast
        iv
    end
    
    methods
        function obj = ivp_var(ast)
            obj.ast = ast;
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for o=obj
                vec = [vec; o.ast.m_value];
            end
        end
        
        function vec = vec(obj)
            vec = [];
            for o=obj
                vec = [vec; o.ast];
            end
        end
        
        function ID = ids(obj)
            ID = [];
            for o = obj
                ID = [ID, o.ast.m_id];
            end
        end
    end
end