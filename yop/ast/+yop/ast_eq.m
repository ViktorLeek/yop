classdef ast_eq < yop.ast_relation
    
    properties
         m_alg = false
    end
    
    properties (Constant)
        name = 'eq'
    end
    
    methods
        function obj = ast_eq(lhs, rhs, ishard, isalg)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = size(eq(ones(size(lhs)), ones(size(rhs))));
            if nargin > 2
                obj.m_hard = ishard;
                obj.m_alg = isalg;
            end
        end
        
        function value = evaluate(obj)
            value = eq(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = eq(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
        
        function obj = alg(obj)
            obj.m_alg = true;
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_eq(lhs, rhs, obj.m_hard, obj.m_alg);
        end
    end
    
end