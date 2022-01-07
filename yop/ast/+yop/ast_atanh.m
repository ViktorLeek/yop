classdef ast_atanh < yop.ast_unary_expression
    properties (Constant)
        name = 'atanh'
    end
    methods
        function obj = ast_atanh(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                atanh(expr.m_value) , ... value
                atanh(expr.m_numval), ... numval
                expr.m_t0           , ... t0
                expr.m_tf           , ... tf
                false(sz)           , ... isder
                expr.m_reducible    , ... isreducible
                zeros(sz)           , ... type
                zeros(sz)           , ... typeid
                expr                 ... expr
                )
        end
    end
end