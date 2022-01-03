classdef ast_independent < yop.ast_variable
    methods
        function obj = ast_independent(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function [bool, id] = isa_independent(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
        end
    end
end