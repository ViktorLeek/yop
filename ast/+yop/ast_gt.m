classdef ast_gt < yop.ast_relation
    
    properties (Constant)
        name = 'gt'
    end
    
    methods
        function obj = ast_gt(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = gt(ones(size(lhs)), ones(size(rhs)));
        end
        
        function value = evaluate(obj)
            value = gt(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = gt(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
    end
    
end