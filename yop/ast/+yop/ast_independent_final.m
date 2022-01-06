classdef ast_independent_final < yop.ast_variable
    methods
        function obj = ast_independent_final(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function [type, id] = Type(obj)
            type = yop.var_type.timef*ones(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end