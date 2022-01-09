classdef ast_ne < yop.ast_relation
    
    properties (Constant)
        m_name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs, ishard)
            if isa(lhs, 'function_handle')
                val = ne(lhs(1), rhs);
            elseif isa(rhs, 'function_handle')
                val = ne(lhs, rhs(1));
            else
                val = ne(value(lhs), value(rhs));
            end
            obj@yop.ast_relation(val, lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_ne(lhs, rhs, obj.m_hard);
        end
    end
    
end