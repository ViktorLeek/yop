classdef ast_mrdivide < yop.ast_binary_expression
    
    methods
        function obj = ast_mrdivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mrdivide(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'mrdivide(lhs, rhs)');
        end
    end
end