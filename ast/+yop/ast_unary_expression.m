classdef ast_unary_expression < yop.ast_expression
    
    properties
        expr
    end
    
    methods
        function obj = ast_unary_expression(expr)
            obj.expr = expr;
        end
        
        function ast(obj)
            fprintf([obj.name, '(expr)\n']);
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
        
    end
end