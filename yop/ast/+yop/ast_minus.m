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
        
        function v = forward(obj)
            obj.m_value = minus(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end