classdef ast_independent_final < yop.ast_variable
    methods
        function obj = ast_independent_final(name)
            obj@yop.ast_variable(name, 1, 1);
        end
        
        function [bool, id] = isa_independent(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end