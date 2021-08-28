classdef ast_independent_initial < yop.ast_variable
    methods
        function obj = ast_independent_initial(name)
            obj@yop.ast_variable(name, 1, 1);
        end
    end
end