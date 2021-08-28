classdef ocp_ode < yop.node
    properties
        var
        ode
    end
    methods
        function obj = ocp_ode(var, ode)
            obj@yop.node();
            obj.var = var;
            obj.ode = ode;
        end
        
        function draw(obj)
            fprintf('ocp_ode(var, ode)\n');
            
            begin_child(obj);
            draw(obj.var);
            end_child(obj);
            
            last_child(obj);
            draw(obj.ode);
            end_child(obj);
        end
    end
end