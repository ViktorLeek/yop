classdef ast_transpose < yop.ast_unary_expression
    properties (Constant)
        m_name = 'transpose'
    end
    methods
        function obj = ast_transpose(expr)
            obj@yop.ast_unary_expression( ...
                transpose(expr.m_value)    , ... value
                transpose(expr.m_numval)   , ... numval
                transpose(expr.m_t0)       , ... t0
                transpose(expr.m_tf)       , ... tf
                transpose(expr.m_der)      , ... der
                transpose(expr.m_reducible), ... reducible
                transpose(expr.m_type)     , ... type
                transpose(expr.m_typeid)   , ... typeid
                expr                         ...
                );
        end
    end
end