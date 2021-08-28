classdef ast_control < yop.ast_variable
    methods
        function obj = ast_control(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
    end
end