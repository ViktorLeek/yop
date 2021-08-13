classdef ast_times < yop.ast_binary_expression
    
    properties (Constant)
        name = 'times'
    end
    
    methods
        function obj = ast_times(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( times(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = times(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end