classdef ast_control < yop.ast_variable
    
    properties
        m_deg
        m_du
    end
    
    methods
        function obj = ast_control(name, w, os, deg)
            obj@yop.ast_variable(name, w, os, false, ...
                (deg>0)*yop.var_type.state + (deg<=0)*yop.var_type.control)
            obj.m_deg = deg;
            if deg > 0
                obj.m_du = yop.ast_control(['D', name], 1, 0, deg-1);
            end
        end
    end
end