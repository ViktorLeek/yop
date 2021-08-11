classdef ast_power < yop.ast_binary_expression
    
    methods
        function obj = ast_power(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( power(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'power(lhs, rhs)');
        end
    end
end