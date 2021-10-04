classdef ast_parameter < yop.ast_variable
    methods
        function obj = ast_parameter(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
        
        function [bool, id] = isa_state(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end