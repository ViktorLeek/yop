classdef ast_ctranspose < yop.ast_unary_expression
    properties (Constant)
        m_name = 'ctranspose'
    end
    methods
        function obj = ast_ctranspose(expr)
            obj@yop.ast_unary_expression( ...
                ctranspose(expr.m_value)    , ... value
                ctranspose(expr.m_numval)   , ... numval
                ctranspose(expr.m_t0)       , ... t0
                ctranspose(expr.m_tf)       , ... tf
                ctranspose(expr.m_der)      , ... der
                ctranspose(expr.m_reducible), ... reducible
                ctranspose(expr.m_type)     , ... type
                ctranspose(expr.m_typeid)   , ... typeid
                expr                         ...
                );
        end
    end
end