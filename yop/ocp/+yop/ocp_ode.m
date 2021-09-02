classdef ocp_ode < handle
    properties
        var
        expr
    end
    methods
        function obj = ocp_ode(var, expr)
            obj.var = var;
            obj.expr = expr;
        end
    end
end