classdef ast_le < yop.ast_relation
    
    properties (Constant)
        m_name = 'le'
    end
    
    methods
        function obj = ast_le(lhs, rhs, ishard)
            if isa(lhs, 'function_handle')
                val = le(lhs(1), value(rhs));
                num = le(lhs(1), numval(rhs));
            elseif isa(rhs, 'function_handle')
                val = le(value(lhs), rhs(1));
                num = le(numval(lhs), rhs(1));
            else
                val = le(value(lhs), value(rhs));
                num = le(numval(lhs), numval(rhs));
            end
            obj@yop.ast_relation(val, num, lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_le(lhs, rhs, obj.m_hard);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.m_lhs-obj.m_rhs, 0);
        end

    end
    
end