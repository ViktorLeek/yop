classdef PI < handle
    % Continuous time PI controller with back calculation for anti-windup
    properties
        K
        Ti
        Tt
    end
    methods
        function obj = PI(K, Ti, Tt)
            obj.K = K;
            obj.Ti = Ti;
            obj.Tt = Tt;
        end
        
        function ctrl = u(obj, e, I)
            ctrl = obj.K*e + I;
        end
        
        function der = dI(obj, e, es)
            % e  = y_sp - y
            % es = u_act - u, u_act may be saturated
            der = obj.K/obj.Ti*e + es/obj.Tt;
        end
    end
end