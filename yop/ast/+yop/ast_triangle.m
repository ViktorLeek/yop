classdef ast_triangle < yop.ast_unary_expression
    properties (Constant)
        m_name = 'triangle'
    end
    methods
        function obj = ast_triangle(expr)
            val = triangle(expr.m_value);
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