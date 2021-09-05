classdef ocp_expr < handle
    properties
        expr
        fn
        disc
    end
    methods
        function obj = ocp_expr(expr)
            obj.expr = expr;
        end
        function bool = is_transcription_invariant(obj)
            bool = is_transcription_invariant(obj.expr);
        end
        function vec = vertcat_disc(obj)
            vec = [];
            for k=1:length(obj)
                vec = [vec(:); obj(k).disc(:)];
            end
        end
    end
end