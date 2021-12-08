classdef ast_control < yop.ast_variable
    properties
        deg
        der
    end
    methods
        function obj = ast_control(name, deg)
            obj@yop.ast_variable(name);
            obj.deg = deg;
            if deg > 0
                obj.der = yop.ast_control(['D', obj.name], deg-1);
            end
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function [bool, id] = isa_control(obj)
            if isempty(obj.der)
                bool = true(size(obj));
            else
                bool = false(size(obj));
            end
            id = obj.id*ones(size(obj));
        end
        
        function [bool, id] = isa_state(obj)
            if isempty(obj.der)
                bool = false(size(obj));
            else
                bool = true(size(obj));
            end
            id = obj.id*ones(size(obj));
        end
    end
end