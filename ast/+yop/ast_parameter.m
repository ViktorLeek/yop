classdef ast_parameter < yop.ast_variable
    methods
        function obj = ast_parameter(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
    end
end