classdef ast_algebraic < yop.ast_variable
    methods
        function obj = ast_algebraic(name, w, os)
            obj@yop.ast_variable(name, w, os, false, yop.var_type.algebraic);
        end
    end
end