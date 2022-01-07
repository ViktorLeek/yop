classdef ast_ceil < yop.ast_unary_expression
    properties (Constant)
        name = 'ceil'
    end
    methods
        function obj = ast_ceil(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                ceil(expr.m_value) , ... value
                ceil(expr.m_numval), ... numval
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