classdef ocp_var < handle
    properties
        var
        mx
        sym
        m_ub0  % initial upper bound
        m_lb0  % initial lower bound
        m_ub   % upper bound
        m_lb   % lower bound
        m_ubf  % final upper bound
        m_lbf  % final lower bound
    end
    properties (Hidden)
        enum
    end
    methods
        function obj = ocp_var(var)       
            obj.var = var;
            obj.mx = yop.cx(['ocp_', var.name], size(var,1), size(var,2));
            obj.sym = sym(var.name, size(var));
            obj.m_ub0 = cell(size(var));
            obj.m_lb0 = cell(size(var));
            obj.m_ub  = cell(size(var));
            obj.m_lb  = cell(size(var));
            obj.m_ubf = cell(size(var));
            obj.m_lbf = cell(size(var));
        end
        
        function v = get_defval(obj, len)
            v = cell(len,1);
            for k=1:len
                v{k} = nan;
            end
        end
        
        function bd = ub0(obj, t)
            bd = get_bnd(obj, 'm_ub0', t);
        end
        
        function bd = lb0(obj, t)
            bd = get_bnd(obj, 'm_lb0', t);
        end
        
        function bd = ub(obj, t)
            bd = get_bnd(obj, 'm_ub', t);
        end
        
        function bd = lb(obj, t)
            bd = get_bnd(obj, 'm_lb', t);
        end
        
        function bd = ubf(obj, t)
            bd = get_bnd(obj, 'm_ubf', t);
        end
        
        function bd = lbf(obj, t)
            bd = get_bnd(obj, 'm_lbf', t);
        end
        
        function bd = get_bnd(obj, prop, t)
            bd = [];
            for o=obj
                for k=1:length(o.(prop))
                    bd = [bd; o.(prop){k}(t)]; % Scalar!
                end
            end
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
        
        function vec = vec(obj)
            vec = [];
            for o=obj
                vec = [vec(:); o.var(:)];
            end
        end
        
        function obj = set_boxfn(obj)
            for o=obj
                for field={'m_ub0','m_lb0','m_ub','m_lb','m_ubf','m_lbf'}
                    for k=1:length(o.(field{1}))
                        if isnumeric(o.(field{1}){k})
                            val = o.(field{1}){k};
                            o.(field{1}){k} = @(t) val;
                        end
                    end
                end
            end
        end
        
    end
end