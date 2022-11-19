classdef ast_state_der < yop.ast_variable
    methods
        function obj = ast_state_der(name, w, os)
            obj@yop.ast_variable(name, w, os, false, yop.var_type.state_der);
        end
    end
end