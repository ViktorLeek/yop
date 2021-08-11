classdef ast_plus < yop.ast_binary_expression
    
    methods
        function obj = ast_plus(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( plus(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'plus(lhs, rhs)');
        end
    end
end