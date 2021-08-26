classdef ast_le < yop.ast_relation
    
    properties (Constant)
        name = 'le'
    end
    
    methods
        function obj = ast_le(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = size(le(ones(size(lhs)), ones(size(rhs))));
        end
        
        function value = evaluate(obj)
            value = le(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = le(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
end