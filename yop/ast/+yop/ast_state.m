classdef ast_state < yop.ast_variable
    methods
        function obj = ast_state(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function [bool, id] = isa_state(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
        
    end
end