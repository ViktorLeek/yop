classdef ocp_independent_final < yop.node
    properties
        var
        ub   % upper bound
        lb   % lower bound
        
    end
    methods
        function obj = ocp_independent_final(var)
            obj@yop.node();
            obj.var = var;
            obj.ub  = yop.defaults().independent_ubf;
            obj.lb  = yop.defaults().independent_lbf;
        end
        
        function draw(obj)
            fprintf('ocp_independent_final(var, ub, lb)\n');
            
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