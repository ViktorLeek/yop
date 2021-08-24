classdef ast_lt < yop.ast_relation
    
    properties (Constant)
        name = 'lt'
    end
    
    methods
        function obj = ast_lt(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = lt(ones(size(lhs)), ones(size(rhs)));
        end
        
        function value = evaluate(obj)
            value = lt(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = lt(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end