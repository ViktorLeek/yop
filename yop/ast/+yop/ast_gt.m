classdef ast_gt < yop.ast_relation
    
    properties (Constant)
        m_name = 'gt'
    end
    
    methods
        function obj = ast_gt(lhs, rhs, ishard)
            if isa(lhs, 'function_handle')
                val = gt(lhs(1), value(rhs));
                num = gt(lhs(1), numval(rhs));
            elseif isa(rhs, 'function_handle')
                val = gt(value(lhs), rhs(1));
                num = gt(numval(lhs), rhs(1));
            else
                val = gt(value(lhs), value(rhs));
                num = gt(numval(lhs), numval(rhs));
            end
            obj@yop.ast_relation(val, num, lhs, rhs);
            if nargin == 3
                obj.m_hard = ishard;
            end
        end
        
        function fn = get_constructor(obj)
            fn = @(lhs, rhs) yop.ast_gt(lhs, rhs, obj.m_hard);
        end
        
        function ceq = canonicalize(obj)
            fn = get_constructor(obj);
            ceq = fn(obj.m_rhs-obj.m_lhs, 0);
        end

    end
    
end