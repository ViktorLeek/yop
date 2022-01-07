classdef ast_atan < yop.ast_unary_expression
    properties (Constant)
        name = 'atan'
    end
    methods
        function obj = ast_atan(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                atan(expr.m_value) , ... value
                atan(expr.m_numval), ... numval
                expr.m_t0          , ... t0
                expr.m_tf          , ... tf
                false(sz)          , ... isder
                expr.m_reducible   , ... isreducible
                zeros(sz)          , ... type
                zeros(sz)          , ... typeid
                expr                ... expr
                )
        end
    end
end