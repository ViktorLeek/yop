classdef ast_asin < yop.ast_unary_expression
    properties (Constant)
        name = 'asin'
    end
    methods
        function obj = ast_asin(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                asin(expr.m_value) , ... value
                asin(expr.m_numval), ... numval
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