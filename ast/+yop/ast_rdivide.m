classdef ast_rdivide < yop.ast_binary_expression
    
    properties (Constant)
        name = 'rdivide'
    end
    
    methods
        function obj = ast_rdivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( rdivide(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = rdivide(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = rdivide(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
end