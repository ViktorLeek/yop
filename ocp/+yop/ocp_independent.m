classdef ocp_independent < yop.node
    properties
        var
        ub   % upper bound
        lb   % lower bound
        
    end
    methods
        function obj = ocp_independent(var)
            obj@yop.node();
            obj.var = var;
            obj.ub  = yop.defaults().independent_ub;
            obj.lb  = yop.defaults().independent_lb;
        end
        
        function draw(obj)
            fprintf('ocp_independent(var, ub, lb)\n');
            
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