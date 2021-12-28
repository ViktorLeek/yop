classdef ast_parameter < yop.ast_variable
    methods
        function obj = ast_parameter(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function [bool, id] = isa_parameter(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end