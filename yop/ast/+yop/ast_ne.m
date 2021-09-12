classdef ast_ne < yop.ast_relation
    
    properties (Constant)
        name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs, ishard)
            obj@yop.ast_relation(lhs, rhs);
            obj.dim = size(ne(ones(size(lhs)), ones(size(rhs))));
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function value = evaluate(obj)
            value = ne(evaluate(obj.lhs), evaluate(obj.rhs));
        end
        
        function v = forward(obj)
            obj.m_value = ne(value(obj.lhs), value(obj.rhs));
            v = obj.m_value;
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_ne(lhs, rhs, obj.m_hard);
        end
    end
    
end