classdef ast_mtimes < yop.ast_binary_expression
    
    properties (Constant)
        name = 'mtimes'
    end
    
    methods
        function obj = ast_mtimes(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mtimes(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = mtimes(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = mtimes(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end