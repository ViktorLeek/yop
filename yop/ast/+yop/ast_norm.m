classdef ast_norm < yop.ast_binary_expression
    
    properties (Constant)
        m_name = 'norm'
    end
    
    methods
        function obj = ast_norm(expr, p)
            if nargin==1
                p = 2;
            end
            num = norm(expr.m_numval, p);
            sz = size(num);
            reducible = all(expr.m_reducible) & true(sz);
            t0 = max(expr.m_t0(:)) * ones(sz);
            tf = min(expr.m_tf(:)) * ones(sz);
            obj@yop.ast_binary_expression( ...
                norm(value(expr), value(p)), ... value
                num                        , ... numval
                t0                         , ... t0
                tf                         , ... tf
                false(sz)                  , ... der
                reducible                  , ... reducible
                zeros(sz)                  , ... type
                zeros(sz)                  , ... typeid
                expr                       , ... lhs
                p                           ... rhs
                );
        end
    end
    
end