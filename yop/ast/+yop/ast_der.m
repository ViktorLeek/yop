classdef ast_der < yop.ast_unary_expression
    
    properties (Constant)
        m_name = 'der'
    end
    
    methods
        function obj = ast_der(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                yop.cx('der', sz(1), sz(2)), ... value
                nan(sz)                    , ... numval
                expr.m_t0                  , ... t0
                expr.m_tf                  , ... tf
                []                         , ... isder
                false(sz)                  , ... reducible
                expr.m_type                , ... type
                expr.m_typeid              , ... typeid
                expr                        ... expr
                );
            obj.m_der = obj.m_id*ones(size(obj)); % A bit ugly
        end
    end
end