classdef ast_ge < yop.ast_relation
    
    properties (Constant)
        name = 'ge'
    end
    
    methods
        function obj = ast_ge(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
        end
        
        function value = evaluate(obj)
            value = ge(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = ge(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end