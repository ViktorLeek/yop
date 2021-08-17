classdef ast_uplus < yop.ast_unary_expression
    
    properties (Constant)
        name = 'uplus'
    end
    
    methods
        
        function obj = ast_uplus(expr)
            obj@yop.ast_unary_expression(expr);
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = uplus(evaluate(obj.expr));
        end
        
        function v = forward(obj)
            obj.m_value = uplus(value(obj.expr));
            v = obj.m_value;
        end
        
    end
    
end