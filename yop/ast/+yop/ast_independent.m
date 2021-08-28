classdef ast_independent < yop.ast_variable
    methods
        function obj = ast_independent(name)
            obj@yop.ast_variable(name, 1, 1);
        end
    end
end