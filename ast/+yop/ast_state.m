classdef ast_state < yop.ast_variable
    methods
        function obj = ast_state(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
    end
end