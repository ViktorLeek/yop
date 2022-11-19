classdef ast_transpose < yop.ast_unary_expression
    properties (Constant)
        m_name = 'transpose'
    end
    properties
        m_dval
    end
    methods
        function obj = ast_transpose(expr)
            obj@yop.ast_unary_expression( ...
                transpose(expr.m_value)    , ... value
                transpose(expr.m_numval)   , ... numval
                transpose(expr.m_t0)       , ... t0
                transpose(expr.m_tf)       , ... tf
                transpose(expr.m_der)      , ... der
                transpose(expr.m_reducible), ... reducible
                transpose(expr.m_type)     , ... type
                transpose(expr.m_typeid)   , ... typeid
                expr                         ...
                );
            if all(obj.m_type == yop.var_type.state)
                obj.m_dval = transpose(expr.m_dval);
            end
        end
        
        function d = der(obj)
            if isempty(obj.m_dval)
                d = der@yop.ast_expression(obj);
            else
                d = obj.m_dval;
            end
        end
    end
end