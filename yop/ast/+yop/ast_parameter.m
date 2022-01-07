classdef ast_parameter < yop.ast_variable
    methods
        function obj = ast_parameter(name, w, os)
            obj@yop.ast_variable(name, w, os, true, yop.var_type.parameter);
        end
    end
end