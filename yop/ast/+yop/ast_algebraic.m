classdef ast_algebraic < yop.ast_variable
    methods
        function obj = ast_algebraic(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function [bool, id] = isa_algebraic(obj)
            bool = true;
            id = obj.id;
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
    end
end