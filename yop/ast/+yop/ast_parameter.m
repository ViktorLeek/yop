classdef ast_parameter < yop.ast_variable
    methods
        function obj = ast_parameter(name)
            obj@yop.ast_variable(name);
        end
        
        function [bool, id] = isa_parameter(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end