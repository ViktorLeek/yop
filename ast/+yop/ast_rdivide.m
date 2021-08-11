classdef ast_rdivide < yop.ast_binary_expression
    
    methods
        function obj = ast_rdivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( rdivide(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'rdivide(lhs, rhs)');
        end
    end
end