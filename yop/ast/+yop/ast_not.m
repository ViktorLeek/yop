classdef ast_not < yop.ast_logical
    
    properties (Constant)
        m_name = 'not'
    end
    
    methods
        function obj = ast_not(lhs, rhs)
            val = not(value(lhs), value(rhs));
            obj@yop.ast_logical(val, lhs, rhs);
        end
    end
    
end