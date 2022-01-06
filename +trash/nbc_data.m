classdef nbc_data < handle
    % Non-box constraint data
    properties
       odes;      % ode constraints  der(x) == f(t,x,z,u,p)
       eq = {};  % equality constraints, g(t,x,z,u,p) == 0
       ieq = {}; % inequality constraints, h(t,x,z,u,p) <= 0 
    end
    methods
        function obj = add_ode(obj, der, expr)
            if isempty(obj.odes)
                obj.odes = yop.ocp_ode(der, expr);
            else
                obj.odes(end+1) = yop.ocp_ode(der, expr);
            end            
        end
        
        function obj = add_ieq(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.ieq = {obj.ieq{:}, e{k}};
                end
            else
                obj.ieq = {obj.ieq{:}, e};
            end
        end
        
        function obj = add_eq(obj, e)
            if isempty(e)
                return
            end
            if iscell(e)
                for k=1:length(e)
                    obj.eq = {obj.eq{:}, e{k}};
                end
            else
                obj.eq = {obj.eq{:}, e};
            end
        end
        
    end
end