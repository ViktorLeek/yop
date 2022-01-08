classdef ast_erf < yop.ast_unary_expression
    properties (Constant)
        m_name = 'erf'
    end
    methods
        function obj = ast_erf(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                erf(expr.m_value) , ... value
                erf(expr.m_numval), ... numval
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