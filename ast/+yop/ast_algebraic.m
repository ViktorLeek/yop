classdef ast_algebraic < yop.ast_variable
    methods
        function obj = ast_algebraic(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
    end
end