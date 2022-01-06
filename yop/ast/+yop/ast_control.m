classdef ast_control < yop.ast_variable
    properties
        deg
        der
    end
    methods
        function obj = ast_control(name, w, os, deg)
            obj@yop.ast_variable(name, w, os);
            obj.deg = deg;
            if deg > 0
                obj.der = yop.ast_control(['D', obj.name], 1, 0, deg-1);
            end
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function [bool, id, type] = isa_variable(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
            if isempty(obj.der)
                type = yop.var_type.control*ones(size(obj));
            else
                type = yop.var_type.state*ones(size(obj));
            end
        end
    end
end