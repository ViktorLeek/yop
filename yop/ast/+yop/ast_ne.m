classdef ast_ne < yop.ast_relation
    
    properties (Constant)
        m_name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs, ishard)
            obj@yop.ast_relation(ne(value(lhs), value(rhs)), lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_ne(lhs, rhs, obj.m_hard);
        end
    end
    
end