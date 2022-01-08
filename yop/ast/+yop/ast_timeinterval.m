classdef ast_timeinterval < yop.ast_unary_expression
    properties (Constant)
        m_name = 'timeinterval'
    end
    methods
        function obj = ast_timeinterval(t0, tf, expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                expr.m_value    , ... value
                expr.m_numval   , ... numval
                t0*ones(sz)     , ... t0
                tf*ones(sz)     , ... tf
                false(sz)       , ... isder
                expr.m_reducible, ... isreducible
                expr.m_type     , ... type
                expr.m_typeid   , ... typeid
                expr             ... expr
                )
        end
    end
end