classdef ast_minus < yop.ast_binary_expression
    
    properties (Constant)
        name = 'minus'
    end
    
    methods
        function obj = ast_minus(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( minus(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = minus(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end