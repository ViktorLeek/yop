classdef ast_or < yop.ast_logical
    
    properties (Constant)
        m_name = 'or'
    end
    
    methods
        function obj = ast_or(lhs, rhs)
            val = or(value(lhs), value(rhs));
            obj@yop.ast_logical(val, lhs, rhs);
        end
    end
    
end