classdef ast_gt < yop.ast_relation
    
    properties (Constant)
        m_name = 'gt'
    end
    
    methods
        function obj = ast_gt(lhs, rhs, ishard)
            obj@yop.ast_relation(gt(value(lhs), value(rhs)), lhs, rhs);
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