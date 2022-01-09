classdef ast_le < yop.ast_relation
    
    properties (Constant)
        m_name = 'le'
    end
    
    methods
        function obj = ast_le(lhs, rhs, ishard)
            %             obj@yop.ast_relation(le(value(lhs), value(rhs)), lhs, rhs);
            obj@yop.ast_relation([], lhs, rhs);
            if nargin > 2
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