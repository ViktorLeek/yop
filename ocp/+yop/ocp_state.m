classdef ocp_state < yop.node
    properties
        var
        ub0  % initial upper bound
        lb0  % initial lower bound
        ub   % upper bound
        lb   % lower bound
        ubf  % final upper bound
        lbf  % final lower bound
        
    end
    methods
        function obj = ocp_state(var)
            obj@yop.node();
            obj.var = var;
            obj.ub0 = yop.defaults().state_ub0 * ones(size(var));
            obj.lb0 = yop.defaults().state_lb0 * ones(size(var));
            obj.ub  = yop.defaults().state_ub  * ones(size(var));
            obj.lb  = yop.defaults().state_lb  * ones(size(var));
            obj.ubf = yop.defaults().state_ubf * ones(size(var));
            obj.lbf = yop.defaults().state_lbf * ones(size(var));
        end
        
        function draw(obj)
            fprintf('ocp_state(var, ub0, lb0, ub, lb, ubf, lbf)\n');
            
            begin_child(obj);
            draw(obj.var);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.ub0);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.lb0);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.ub);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.lb);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.ubf);
            end_child(obj);
            
            last_child(obj);
            draw(obj.lbf);
            end_child(obj);
        end
    end
end