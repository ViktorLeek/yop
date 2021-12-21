classdef ivp_var < handle
    properties
        var
        mx
        sym
        ub
        lb
    end
    properties (Hidden)
        enum
    end
    methods
        function obj = ivp_var(var)       
            obj.var = var;
            obj.mx = yop.cx(['ivp_', var.name], size(var,1), size(var,2));
            obj.sym = sym(var.name, size(var));
            obj.ub = nan(size(var));
            obj.lb = nan(size(var));
        end
        
        function n = n_elem(obj)
            n = 0;
            for o=obj
                n = n + prod(size(o.var));
            end
        end
        
        function obj = set_value(obj, value)
            obj.var.m_value = value;
        end
        
        function obj = set_sym(obj)
            for o=obj
                o.var.m_value = o.sym;
            end
        end
        
        function obj = set_mx(obj)
            for o=obj
                o.var.m_value = o.mx;
            end
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for o=obj
                vec = [vec(:); o.mx(:)];
            end
        end
        
        function ids = IDs(obj)
            ids = [];
            for o=obj
                ids = [ids, o.var.id];
            end
        end
    end
end