classdef ast_plus < yop.ast_binary_expression
    
    properties (Constant)
        m_name = 'plus'
    end
    
    methods
        function obj = ast_plus(lhs, rhs)
            num = plus(numval(lhs), numval(rhs));
            sz = size(num);
            obj@yop.ast_binary_expression( ...
                plus(value(lhs), value(rhs)) , ... value
                num                          , ... numval
                max(get_t0(lhs), get_t0(rhs)), ... t0
                min(get_tf(lhs), get_tf(lhs)), ... tf
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