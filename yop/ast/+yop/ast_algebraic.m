classdef ast_algebraic < yop.ast_variable
    methods
        function obj = ast_algebraic(name)
            obj@yop.ast_variable(name);
        end
        
        function bool = is_algebraic(obj)
            bool = true;
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
    end
end