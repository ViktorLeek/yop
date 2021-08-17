classdef ast_independent_final < yop.ast_variable
    methods
        function obj = ast_independent_final(name)
            obj@yop.ast_variable(name, 1, 1);
        end
    end
end