classdef ast_ge < yop.ast_relation
    
    properties (Constant)
        m_name = 'ge'
    end
    
    methods
        function obj = ast_ge(lhs, rhs, ishard)
            if isa(lhs, 'function_handle')
                val = ge(lhs(1), rhs);
            elseif isa(rhs, 'function_handle')
                val = ge(lhs, rhs(1));
            else
                val = ge(value(lhs), value(rhs));
            end
            obj@yop.ast_relation(val, lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_ge(lhs, rhs, obj.m_hard);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.m_rhs-obj.m_lhs, 0);
        end

    end
    
end