classdef ast_control < yop.ast_variable
    methods
        function obj = ast_control(name, rows, cols)
            obj@yop.ast_variable(name, rows, cols);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
    end
end