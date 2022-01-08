classdef ast_power < yop.ast_binary_expression
    
    properties (Constant)
        m_name = 'power'
    end
    
    methods
        function obj = ast_power(lhs, rhs)
            num = power(numval(lhs), numval(rhs));
            sz = size(num);
            t0_l = get_t0(lhs);
            t0_r = get_t0(rhs);
            tf_l = get_tf(lhs);
            tf_r = get_tf(rhs);
            t0 = max([t0_l(:); t0_r(:)]) * ones(sz);
            tf = min([tf_l(:); tf_r(:)]) * ones(sz);
            obj@yop.ast_binary_expression( ...
                power(value(lhs), value(rhs)) , ... value
                num                          , ... numval
                t0                           , ... t0
                tf                           , ... tf
                false(sz)                    , ... der
                isa_reducible(lhs) & isa_reducible(rhs) , ... reducible
                zeros(sz)                    , ... type
                zeros(sz)                    , ... typeid
                lhs                          , ... lhs
                rhs                           ... rhs
                );
        end
    end
    
end