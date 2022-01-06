classdef ast_independent_initial < yop.ast_variable
    methods
        function obj = ast_independent_initial(name, w, os)
            obj@yop.ast_variable(name, w, os);
            obj.m_value = yop.cx(name);
        end
        
        function [type, id] = Type(obj)
            type = yop.var_type.time0*ones(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end