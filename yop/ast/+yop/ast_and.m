classdef ast_and < yop.ast_logical
    
    properties (Constant)
        m_name = 'and'
    end
    
    methods
        function obj = ast_and(lhs, rhs)
            val = and(value(lhs), value(rhs));
            obj@yop.ast_logical(val, lhs, rhs);
        end
    end
    
end