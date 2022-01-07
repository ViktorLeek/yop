classdef ast_der < yop.ast_unary_expression
    
    properties (Constant)
        name = 'der'
    end
    
    methods
        function obj = ast_der(expr)
            sz = size(expr);
            obj@yop.ast_unary_expression( ...
                yop.cx('der', sz(1), sz(2)), ... value
                nan(sz)                    , ... numval
                yop.initial_timepoint(sz)  , ... t0
                yop.final_timepoint(sz)    , ... tf
                []                         , ... isder
                false(sz)                  , ... reducible
                zeros(sz)                  , ... type
                zeros(sz)                  , ... typeid
                expr                        ... expr
                );
            obj.m_der = obj.id*ones(size(obj)); % ugly fix
        end
    end
end