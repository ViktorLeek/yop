classdef ast_mldivide < yop.ast_binary_expression
    
    properties (Constant)
        name = 'mldivide'
    end
    
    methods
        function obj = ast_mldivide(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mldivide(ones(size(lhs)), ones(size(rhs))) );
        end
        
        function value = evaluate(obj)
            value = mldivide(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = mldivide(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
   
end