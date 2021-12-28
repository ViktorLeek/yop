classdef ast_independent_initial < yop.ast_variable
    methods
        function obj = ast_independent_initial(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function [bool, id] = isa_independent0(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end