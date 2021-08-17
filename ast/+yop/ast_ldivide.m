classdef ast_ldivide < yop.ast_binary_expression
    
    properties (Constant)
        name = 'ldivide'
    end
    
    methods
        function obj = ast_ldivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( ldivide(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = ldivide(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = ldivide(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
end