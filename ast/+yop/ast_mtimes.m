classdef ast_mtimes < yop.ast_binary_expression
    
    methods
        function obj = ast_mtimes(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mtimes(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'mtimes(lhs, rhs)');
        end
    end
end