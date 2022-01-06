classdef ast_algebraic < yop.ast_variable
    methods
        function obj = ast_algebraic(name, w, os)
            obj@yop.ast_variable(name, w, os);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function [bool, id, type] = isa_variable(obj)
            bool = true(size(obj));
            id = obj.id*ones(size(obj));
            type = yop.var_type.algebraic*ones(size(obj));
        end
    end
end