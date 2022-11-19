classdef ast_control < yop.ast_variable
    
    properties
        m_deg
        m_du
        m_dval
    end
    
    methods
        function obj = ast_control(name, w, os, deg)
            obj@yop.ast_variable(name, w, os, false, ...
                (deg>0)*yop.var_type.state + (deg<=0)*yop.var_type.control)
            obj.m_deg = deg;
            if deg > 0
                obj.m_du = yop.ast_control(['D', name], 1, 0, deg-1);
                obj.m_dval = yop.ast_state_der(['D', name], 1+0*w, 0*os);
            end
        end
        
        function d = du(obj)
            d = obj.m_du;
        end
        
        function du = der(obj)
            if isempty(obj.m_du)
                du = der@yop.ast_variable(obj);
            else
                du = obj.m_du;
            end
        end
    end
end