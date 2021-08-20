classdef ocp_parameter < yop.node
    properties
        var
        ub   % upper bound
        lb   % lower bound
    end
    methods
        function obj = ocp_parameter(var)
            obj@yop.node();
            obj.var = var;
            obj.ub  = yop.defaults().parameter_ub  * ones(size(var));
            obj.lb  = yop.defaults().parameter_lb  * ones(size(var));
        end
        
        function draw(obj)
            fprintf('ocp_parameter(var, ub, lb)\n');
            
            begin_child(obj);
            draw(obj.var);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.ub);
            end_child(obj);
            
            last_child(obj);
            draw(obj.lb);
            end_child(obj);
        end
    end
end