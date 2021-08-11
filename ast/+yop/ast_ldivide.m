classdef ast_ldivide < yop.ast_binary_expression
    
    methods
        function obj = ast_ldivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( ldivide(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'ldivide(lhs, rhs)');
        end
    end
end