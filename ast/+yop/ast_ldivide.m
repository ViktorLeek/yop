classdef ast_ldivide < yop.ast_binary_expression
    
    methods
        function obj = ast_ldivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( ldivide(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = ldivide(evaluate(obj.lhs), evaluate(obj.rhs));
        end

        function ast(obj)
            xast(obj, 'ldivide(lhs, rhs)');
        end
    end
end