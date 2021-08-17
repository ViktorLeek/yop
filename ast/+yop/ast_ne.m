classdef ast_ne < yop.ast_relation
    
    properties (Constant)
        name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
        end
        
        function value = evaluate(obj)
            value = ne(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = ne(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end