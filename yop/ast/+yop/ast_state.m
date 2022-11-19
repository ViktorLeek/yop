classdef ast_state < yop.ast_variable
    properties
        m_dval
    end
    methods
        function obj = ast_state(name, w, os)
            obj@yop.ast_variable(name, w, os, false, yop.var_type.state);
            obj.m_dval = yop.ast_state_der(['D', name], 1+0*w, 0*os);
        end
        
        function dx = der(obj)
            dx = obj.m_dval;
        end
        
        function d = dx(obj)
            d = obj.m_dval;
        end
    end
end