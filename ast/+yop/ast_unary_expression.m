classdef ast_unary_expression < yop.ast_expression
    
    properties
        expr
    end
    
    methods
        function obj = ast_unary_expression(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
        end
        
        function draw(obj)
            fprintf([obj.name, '(expr)\n']);
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
        
    end
end