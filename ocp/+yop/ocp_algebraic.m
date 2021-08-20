classdef ocp_algebraic < yop.node
    properties
        var
        ub   % upper bound
        lb   % lower bound
    end
    methods
        function obj = ocp_algebraic(var)
            obj@yop.node();
            obj.var = var;
            obj.ub  = yop.defaults().algebraic_ub  * ones(size(var));
            obj.lb  = yop.defaults().algebraic_lb  * ones(size(var));
        end
        
        function draw(obj)
            fprintf('ocp_algebraic(var, ub, lb)\n');
            
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