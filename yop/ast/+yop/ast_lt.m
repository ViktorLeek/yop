classdef ast_lt < yop.ast_relation
    
    properties (Constant)
        m_name = 'lt'
    end
    
    methods
        function obj = ast_lt(lhs, rhs, ishard)
            if isa(lhs, 'function_handle')
                val = lt(lhs(1), rhs);
            elseif isa(rhs, 'function_handle')
                val = lt(lhs, rhs(1));
            else
                val = lt(value(lhs), value(rhs));
            end
            obj@yop.ast_relation(val, lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_lt(lhs, rhs, obj.m_hard);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.m_lhs-obj.m_rhs, 0);
        end
        
    end
    
end