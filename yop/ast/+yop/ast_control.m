classdef ast_control < yop.ast_variable
    properties
        du
    end
    methods
        function obj = ast_control(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
        
        function du = get_du(obj)
            if isempty(obj.du)
                obj.du = yop.ast_control(['D', obj.name], size(obj,1), size(obj,2));
            end
            du = obj.du;
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function [bool, id] = isa_control(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end