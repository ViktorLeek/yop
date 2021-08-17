classdef ast_plus < yop.ast_binary_expression
    
    properties (Constant)
        name = 'plus'
    end
    
    methods
        function obj = ast_plus(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( plus(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = plus(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = plus(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end