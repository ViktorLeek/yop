classdef ast_power < yop.ast_binary_expression
    
    properties (Constant)
        name = 'power'
    end
    
    methods
        function obj = ast_power(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( power(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = power(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = power(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end