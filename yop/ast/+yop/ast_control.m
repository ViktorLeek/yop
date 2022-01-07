classdef ast_control < yop.ast_variable
    
    properties (Constant)
        m_reducible = false;
    end
    
    properties
        deg
        der
    end
    
    methods
        function obj = ast_control(name, w, os, deg)
            obj@yop.ast_variable(name, w, os, false, ...
                (deg>0)*yop.var_type.state + (deg<=0)*yop.var_type.control)
            obj.deg = deg;
            if deg > 0
                obj.der = yop.ast_control(['D', obj.name], 1, 0, deg-1);
            end
        end
    end
end