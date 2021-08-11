classdef ast_mldivide < yop.ast_binary_expression
    
    methods
        function obj = ast_mldivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mldivide(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'mldivide(lhs, rhs)');
        end
    end
end