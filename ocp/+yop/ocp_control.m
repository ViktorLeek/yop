classdef ocp_control < yop.node
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
        function obj = ocp_control(var)
            obj@yop.node();
            obj.var = var;
            obj.ub0 = yop.defaults().control_ub0 * ones(size(var));
            obj.lb0 = yop.defaults().control_lb0 * ones(size(var));
            obj.ub  = yop.defaults().control_ub  * ones(size(var));
            obj.lb  = yop.defaults().control_lb  * ones(size(var));
            obj.ubf = yop.defaults().control_ubf * ones(size(var));
            obj.lbf = yop.defaults().control_lbf * ones(size(var));
        end
        
        function draw(obj)
            fprintf('ocp_control(var, ub0, lb0, ub, lb, ubf, lbf)\n');
            
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