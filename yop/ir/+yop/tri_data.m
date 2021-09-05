classdef tri_data < handle
    % trancscription invariant data
    properties
        eq = yop.ocp_expr.empty(1,0)  % equality constraint, transcription invariant
        ieq = yop.ocp_expr.empty(1,0) % inequality constraint, transcription invariant
    end
    methods
        
        function obj = add_eq(obj, e)
            obj.eq(end+1) = yop.ocp_expr(e.lhs);
        end
        
        function obj = add_ieq(obj, e)
            obj.eq(end+1) = yop.ocp_expr(e.lhs);
        end
        
    end
end