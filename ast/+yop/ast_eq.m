classdef ast_eq < yop.ast_relation
    
    properties (Constant)
        name = 'eq'
    end
    
    methods
        function obj = ast_eq(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = size(eq(ones(size(lhs)), ones(size(rhs))));
        end
        
        function value = evaluate(obj)
            value = eq(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = eq(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end