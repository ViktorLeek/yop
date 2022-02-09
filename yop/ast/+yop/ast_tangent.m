classdef ast_tangent < yop.ast_binary_expression
    
    properties (Constant)
        m_name = 'tangent';
    end
    
    methods
        function obj = ast_tangent(lhs, rhs)
            val = tangent(value(lhs), value(rhs));
            sz = size(val);
            num = zeros(sz);
            obj@yop.ast_binary_expression( ...
                val      , ... value
                num      , ... numval
                -inf(sz), ... t0
                +inf(sz), ... tf
                false(sz), ... der
                false(sz), ... reducible
                zeros(sz), ... type
                zeros(sz), ... typeid
                lhs      , ... lhs
                rhs       ... rhs
                );
        end
    end
    
end