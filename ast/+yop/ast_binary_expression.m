classdef ast_binary_expression < yop.ast_node
    
    properties
        lhs
        rhs
    end
    
    methods
        function obj = ast_binary_expression(lhs, rhs)
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
    end
    
    methods % Printing
        
        function xast(obj, s)
            % xast - print with string s as name
            
            
            fprintf([s, '\n']);
            
            begin_child(obj);
            ast(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            ast(obj.rhs);
            end_child(obj);
        end
        
    end
end