classdef ast_sumsqr < yop.ast_unary_expression
    properties (Constant)
        name = 'sumsqr'
    end
    methods
        function obj = ast_sumsqr(expr)
            numval = sumsqr(expr.m_numval);
            sz = size(numval);
            reducible = all(expr.m_reducible) & true(sz);
            t0 = max(expr.m_t0(:)) * ones(sz);
            tf = min(expr.m_tf(:)) * ones(sz);
            obj@yop.ast_unary_expression( ...
                sumsqr(expr.m_value), ... value
                numval           , ... numval
                t0               , ... t0
                tf               , ... tf
                false(sz)        , ... isder
                reducible        , ... isreducible
                zeros(sz)        , ... type
                zeros(sz)        , ... typeid
                expr              ... expr
                )
        end
    end
end