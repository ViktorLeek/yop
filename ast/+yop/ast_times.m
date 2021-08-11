classdef ast_times < yop.ast_binary_expression
    
    methods
        function obj = ast_times(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( times(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'times(lhs, rhs)');
        end
    end
end