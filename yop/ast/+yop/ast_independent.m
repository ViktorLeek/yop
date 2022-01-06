classdef ast_independent < yop.ast_variable
    methods
        function obj = ast_independent(name, w, os)
            obj@yop.ast_variable(name, w, os);
            obj.m_value = yop.cx(name);
        end
        
        function boolv = isa_reducible(obj)
            boolv = false(size(obj));
        end
        
        function [type, id] = Type(obj)
            type = yop.var_type.time*ones(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end