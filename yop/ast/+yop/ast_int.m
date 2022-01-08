classdef ast_int < yop.ast_unary_expression
    
    properties (Constant)
        m_name = 'int'
    end
    
    methods
        function obj = ast_int(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                yop.cx('int', sz(1), sz(2)), ... value
                nan(sz)                    , ... numval
                yop.initial_timepoint(sz)  , ... t0
                yop.final_timepoint(sz)    , ... tf
                false(sz)                  , ... isder
                true(sz)                   , ... reducible
                zeros(sz)                  , ... type
                zeros(sz)                  , ... typeid
                expr                        ...
                );
        end
    end
end