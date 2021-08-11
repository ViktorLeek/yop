classdef ast_unary_expression < yop.ast_node
    
    properties
        expr
    end
    
    methods
        function obj = ast_unary_expression(expr)
            obj.expr = expr;
        end
    end
    
    methods % Printing
        
        function xast(obj, s)
            fprintf([s, '\n']);
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
        
    end
end