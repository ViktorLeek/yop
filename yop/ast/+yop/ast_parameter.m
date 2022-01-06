classdef ast_parameter < yop.ast_variable
    methods
        function obj = ast_parameter(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function [type, id] = Type(obj)
            type = yop.var_type.parameter*ones(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end