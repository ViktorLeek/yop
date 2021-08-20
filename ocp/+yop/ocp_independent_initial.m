classdef ocp_independent_initial < yop.node
    properties
        var
        ub   % upper bound
        lb   % lower bound
        
    end
    methods
        function obj = ocp_independent_initial(var)
            obj@yop.node();
            obj.var = var;
            obj.ub  = yop.defaults().independent_ub0;
            obj.lb  = yop.defaults().independent_lb0;
        end
        
        function draw(obj)
            fprintf('ocp_independent_initial(var, ub, lb)\n');
            
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