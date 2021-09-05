classdef ocp_var < handle
    properties
        var
        mx
        sym
        ub0  % initial upper bound
        lb0  % initial lower bound
        ub   % upper bound
        lb   % lower bound
        ubf  % final upper bound
        lbf  % final lower bound
    end
    properties (Hidden)
        enum
    end
    methods
        function obj = ocp_var(var)       
            obj.var = var;
            obj.mx = casadi.MX.sym(['ocp_', var.name], size(var,1), size(var,2));
            obj.sym = sym(var.name, size(var));
            obj.ub0 = nan(size(var));
            obj.lb0 = nan(size(var));
            obj.ub  = nan(size(var));
            obj.lb  = nan(size(var));
            obj.ubf = nan(size(var));
            obj.lbf = nan(size(var));
        end
        
        function obj = set_value(obj, value)
            obj.var.m_value = value;
        end
        
        
        function obj = set_sym(obj)
            for k=1:length(obj)
                obj(k).var.m_value = obj(k).sym;
            end
        end
        
        function obj = set_mx(obj)
            for k=1:length(obj)
                obj(k).var.m_value = obj(k).mx;
            end
        end
        
        function e = get_enum(obj)
            e = [];
            for k=1:length(obj)
                e = [e(:); obj(k).enum(:)];
            end
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for k=1:length(obj)
                vec = [vec(:); obj(k).mx(:)];
            end
        end
    end
end