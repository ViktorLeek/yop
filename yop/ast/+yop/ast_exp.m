classdef ast_exp < yop.ast_unary_expression
    properties (Constant)
        m_name = 'exp'
    end
    methods
        function obj = ast_exp(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                exp(expr.m_value) , ... value
                exp(expr.m_numval), ... numval
                expr.m_t0         , ... t0
                expr.m_tf         , ... tf
                false(sz)         , ... isder
                expr.m_reducible  , ... isreducible
                zeros(sz)         , ... type
                zeros(sz)         , ... typeid
                expr               ... expr
                )
        end
    end
end