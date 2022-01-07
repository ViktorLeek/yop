classdef ast_state < yop.ast_variable
    methods
        function obj = ast_state(name, w, os)
            obj@yop.ast_variable(name, w, os, false, yop.var_type.state);
        end
    end
end