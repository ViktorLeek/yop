classdef ast_le < yop.ast_relation
    
    properties (Constant)
        name = 'le'
    end
    
    methods
        function obj = ast_le(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
        end
        
        function value = evaluate(obj)
            value = le(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end