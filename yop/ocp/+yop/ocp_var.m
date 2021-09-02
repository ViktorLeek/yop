classdef ocp_var < handle
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
        function obj = ocp_var(var)       
            obj.var = var;
            obj.sym = casadi.MX.sym(var.name, size(var,1), size(var,2));
            obj.ub0 = nan(size(var));
            obj.lb0 = nan(size(var));
            obj.ub  = nan(size(var));
            obj.lb  = nan(size(var));
            obj.ubf = nan(size(var));
            obj.lbf = nan(size(var));
        end
        
        function obj = set_value(obj, value)
            % Push the old value onto the stack
            obj.tmp_value{end+1} = obj.var.m_value;
            obj.var.m_value = value;
        end
        
        function obj = reset_value(obj)
            % Set the value to stack top and pop the stack
            if isempty(obj.tmp_value)
                obj.var.m_value = [];
                obj.tmp_value = {};
            else
                obj.var.m_value = obj.tmp_value{end};
                if length(obj.tmp_value) == 1
                    obj.tmp_value = {};
                else
                    obj.tmp_value = obj.tmp_value{1:end-1};
                end
            end
        end
        
        function obj = set_sym(obj)
            % Push the old value onto the stack
            %obj.tmp_value{end+1} = obj.var.m_value;
            obj.var.m_value = obj.sym;
        end
    end
end