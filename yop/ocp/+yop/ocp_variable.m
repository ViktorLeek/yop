classdef ocp_variable < yop.node
    properties
        var
        sym
        ub0  % initial upper bound
        lb0  % initial lower bound
        ub   % upper bound
        lb   % lower bound
        ubf  % final upper bound
        lbf  % final lower bound
    end
    properties (Hidden)
        tmp_value = {}
    end
    methods
        function obj = ocp_variable(var)
            obj@yop.node();            
            obj.var = var;
            obj.sym = casadi.MX.sym(var.name, size(var,1), size(var,2));
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
            % Push the old value onto the stack
            obj.tmp_value{end+1} = obj.var.m_value;
            obj.var.m_value = value;
        end
        
        function obj = reset_value(obj)
            % Set the value to stack top and pop the stack
            obj.var.m_value = obj.tmp_value{end};
            if length(obj.tmp_value) == 1
                obj.tmp_value = {};
            else
                obj.tmp_value = obj.tmp_value{1:end-1};
            end
        end
        
        function obj = set_sym(obj)
            % Push the old value onto the stack
            obj.tmp_value{end+1} = obj.var.m_value;
            obj.var.m_value = obj.sym;
        end
    end
end