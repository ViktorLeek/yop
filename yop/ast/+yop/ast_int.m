classdef ast_int < yop.ast_unary_expression
    
    properties (Constant)
        m_name = 'int'
    end
    
    properties
        m_to
        m_fr
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
        
        function obj = fr(obj, t)
            obj.m_fr = t;
        end
        
        function obj = to(obj, t)
            obj.m_to = t;
        end
    end
end