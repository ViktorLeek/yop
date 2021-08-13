classdef ast_mpower < yop.ast_binary_expression
    
    properties (Constant)
        name = 'mpower'
    end
    
    methods
        function obj = ast_mpower(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mpower(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = mpower(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
   
end