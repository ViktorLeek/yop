classdef ast_mrdivide < yop.ast_binary_expression
    
    properties (Constant)
        name = 'mrdivide'
    end
    
    methods
        function obj = ast_mrdivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mrdivide(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = mrdivide(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end