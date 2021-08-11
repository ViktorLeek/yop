classdef ast_mpower < yop.ast_binary_expression
    
    methods
        function obj = ast_mpower(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mpower(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'mpower(lhs, rhs)');
        end
    end
end