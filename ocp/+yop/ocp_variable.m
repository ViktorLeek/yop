classdef ocp_variable < yop.node
    properties
        var
        ub0  % initial upper bound
        lb0  % initial lower bound
        ub   % upper bound
        lb   % lower bound
        ubf  % final upper bound
        lbf  % final lower bound
    end
    properties (Hidden)
        tmp_value
    end
    methods
        function obj = ocp_variable(var)
            obj@yop.node();
            obj.var = var;
            obj.ub0 = nan(size(var));
            obj.lb0 = nan(size(var));
            obj.ub  = nan(size(var));
            obj.lb  = nan(size(var));
            obj.ubf = nan(size(var));
            obj.lbf = nan(size(var));
        end
        
        function draw(obj)
            fprintf('ocp_variable(var, ub0, lb0, ub, lb, ubf, lbf)\n');
            
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
        
        function obj = set_value(obj, value)
            obj.var.m_value = value;
        end
        
        function obj = store_value(obj)
            obj.tmp_value = obj.var.m_value;
        end
        
        function obj = restore_value(obj)
            obj.var.m_value = obj.tmp_value;
        end
    end
end