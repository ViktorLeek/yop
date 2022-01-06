classdef ast_algebraic < yop.ast_variable
    methods
        function obj = ast_algebraic(name, w, os)
            obj@yop.ast_variable(name, w, os);
            obj.m_value = yop.cx(name);
        end
        
        function boolv = isa_reducible(obj)
            boolv = false(size(obj));
        end
        
        function [type, id] = Type(obj)
            type = yop.var_type.algebraic*ones(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end