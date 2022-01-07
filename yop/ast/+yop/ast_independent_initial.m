classdef ast_independent_initial < yop.ast_variable
    methods
        function obj = ast_independent_initial(name, w, os)
            obj@yop.ast_variable(name, w, os, true, yop.var_type.time0);
        end
    end
end