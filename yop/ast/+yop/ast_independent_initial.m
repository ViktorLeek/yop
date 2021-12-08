classdef ast_independent_initial < yop.ast_variable
    methods
        function obj = ast_independent_initial(name)
            obj@yop.ast_variable(name);
        end
        
        function [bool, id] = isa_independent(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end