classdef ast_log10 < yop.ast_unary_expression
    properties (Constant)
        m_name = 'log10'
    end
    methods
        function obj = ast_log10(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                log10(expr.m_value) , ... value
                log10(expr.m_numval), ... numval
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