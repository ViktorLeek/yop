classdef ast_independent_final < yop.ast_variable
    methods
        function obj = ast_independent_final(name, w, os)
            obj@yop.ast_variable(name, w, os, true, yop.var_type.timef);
        end
    end
end