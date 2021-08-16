classdef ast_ge < yop.ast_relation
    
    properties (Constant)
        name = 'ge'
    end
    
    methods
        function obj = ast_ge(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
        end
        
        function value = evaluate(obj)
            value = ge(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end