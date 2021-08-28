classdef ast_uminus < yop.ast_unary_expression
    
    properties (Constant)
        name = 'uminus'
    end
    
    methods
        
        function obj = ast_uminus(expr)
            obj@yop.ast_unary_expression(expr);
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = uminus(evaluate(obj.expr));
        end
        
        function v = forward(obj)
            obj.m_value = uminus(value(obj.expr));
            v = obj.m_value;
        end
        
    end
    
end