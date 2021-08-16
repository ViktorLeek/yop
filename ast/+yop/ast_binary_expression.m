classdef ast_binary_expression < yop.ast_expression
    
    properties
        lhs
        rhs
    end
    
    methods
        function obj = ast_binary_expression(lhs, rhs)
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
        
        function ast(obj)
            fprintf([obj.name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            ast(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            ast(obj.rhs);
            end_child(obj);
        end
        
    end
end