classdef ocp_rel < handle
    properties        
        ast
        fn
    end
    methods
        function obj = ocp_rel(rel)
            obj.ast = rel;
        end
        
        function l = lhs(obj)
            l = obj.ast.lhs;
        end
        
        function r = rhs(obj)
            r = obj.ast.rhs;
        end
        
    end
end