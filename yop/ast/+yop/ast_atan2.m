classdef ast_atan2 < yop.ast_binary_expression
    
    properties (Constant)
        m_name = 'atan2'
    end
    
    methods
        function obj = ast_atan2(lhs, rhs)
            num = atan2(numval(lhs), numval(rhs));
            sz = size(num);
            reducible = all(isa_reducible(lhs) & isa_reducible(rhs)) & true(sz);
            obj@yop.ast_binary_expression( ...
                atan2(value(lhs), value(rhs)), ... value
                num                          , ... numval
                max(get_t0(lhs), get_t0(rhs)), ... t0
                min(get_tf(lhs), get_tf(lhs)), ... tf
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