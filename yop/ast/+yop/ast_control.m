classdef ast_control < yop.ast_variable
    
    properties (Constant)
        m_reducible = false;
    end
    
    properties
        m_deg
        m_der
    end
    
    methods
        function obj = ast_control(name, w, os, deg)
            obj@yop.ast_variable(name, w, os, false, ...
                (deg>0)*yop.var_type.state + (deg<=0)*yop.var_type.control)
            obj.m_deg = deg;
            if deg > 0
                obj.m_der = yop.ast_control(['D', name], 1, 0, deg-1);
            end
        end
    end
end