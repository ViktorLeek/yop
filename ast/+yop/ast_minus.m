classdef ast_minus < yop.ast_binary_expression
    
    methods
        function obj = ast_minus(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( minus(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'minus(lhs, rhs)');
        end
    end
end