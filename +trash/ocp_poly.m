classdef ocp_poly < yop.lagrange_polynomial
    properties
        t0
        tf
    end
    
    methods
        function obj = ocp_poly(x, y, t0, tf)
            obj@yop.lagrange_polynomial(x, y);
            obj.t0 = t0;
            obj.tf = tf;
        end
        
        function v = value(obj, x)
            if isa(x, 'yop.ast_independent_initial')
                v = obj(1).evaluate(0);
                return
                
            elseif isa(x, 'yop.ast_independent_final')
                v = obj(end).evaluate(0);
                return
            end
                
            h = (obj.tf - obj.t0)/(length(obj)-1);
            if yop.EQ(rem(x-obj.t0, h), 0)
                n = 1 + round((x-obj.t0)/h);
                v = obj(n).evaluate(0);
                
            else
                n = 1 + floor((x-obj.t0)/h);
                tau = rem(x-obj.t0, h)/h;
                v = obj(n).evaluate(tau);
            end
        end
        
        function y = vec(obj)
            y = [];
            for obj_k = obj(1:end)
                y = [y; obj_k.y(:)];
            end
        end
        
        function y = mat(obj)
            y = [];
            for obj_k = obj(1:end)
                y = [y, obj_k.y];
            end
        end
        
    end
end