classdef ast_rectangle < yop.ast_unary_expression
    properties (Constant)
        m_name = 'rectangle'
    end
    methods
        function obj = ast_rectangle(expr)
            val = rectangle(expr.m_value);
            sz = size(val);
            obj@yop.ast_unary_expression( ...
                val             , ... value
                nan(sz)         , ... numval
                expr.m_t0       , ... t0
                expr.m_tf       , ... tf
                false(sz)       , ... isder
                expr.m_reducible, ... isreducible
                zeros(sz)       , ... type
                zeros(sz)       , ... typeid
                expr             ... expr
                )
        end
    end
end