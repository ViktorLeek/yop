classdef unary_expression < yop.node
    
    properties
        expr
    end
    
    methods
        function obj = binary_expression(expr)
            obj.expr = expr;
        end
    end
    
    methods % Printing
        
        function xprint(obj, s)
            fprintf([s, '\n']);
            last_child(obj);
            yop.print(obj.expr);
            end_child(obj);
        end
        
    end
end