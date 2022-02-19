classdef ast_ne < yop.ast_relation
    
    properties (Constant)
        m_name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs, ishard)
            if isa(lhs, 'function_handle')
                val = ne(lhs(1), value(rhs));
                num = ne(lhs(1), numval(rhs));
            elseif isa(rhs, 'function_handle')
                val = ne(value(lhs), rhs(1));
                num = ne(numval(lhs), rhs(1));
            else
                val = ne(value(lhs), value(rhs));
                num = ne(numval(lhs), numval(rhs));
            end
            obj@yop.ast_relation(val, num, lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_ne(lhs, rhs, obj.m_hard);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.m_rhs-obj.m_lhs, 0);
        end

    end
    
end