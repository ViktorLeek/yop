classdef ast_mpower < yop.ast_binary_expression
    
    properties (Constant)
        m_name = 'mpower'
    end
    
    methods
        function obj = ast_mpower(lhs, rhs)
            num = mpower(numval(lhs), numval(rhs));
            sz = size(num);
            reducible = all(isa_reducible(lhs) & isa_reducible(rhs)) & true(sz);
            t0_l = get_t0(lhs);
            t0_r = get_t0(rhs);
            tf_l = get_tf(lhs);
            tf_r = get_tf(rhs);
            t0 = max([t0_l(:); t0_r(:)]) * ones(sz);
            tf = min([tf_l(:); tf_r(:)]) * ones(sz);
            obj@yop.ast_binary_expression( ...
                mpower(value(lhs), value(rhs)), ... value
                num                          , ... numval
                t0                           , ... t0
                tf                           , ... tf
                false(sz)                    , ... der
                reducible                    , ... reducible
                zeros(sz)                    , ... type
                zeros(sz)                    , ... typeid
                lhs                          , ... lhs
                rhs                           ... rhs
                );
        end
    end
    
end