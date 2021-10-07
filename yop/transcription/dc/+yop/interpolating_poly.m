classdef interpolating_poly < yop.lagrange_polynomial
    properties
        %exclude_last = false;
        N % Number of segments, not necessarily the same as no. elements.
        t0
        tf
    end
    
    methods
        function obj = interpolating_poly(x, y, t0, tf, N)
            obj@yop.lagrange_polynomial(x, y);
            obj.t0 = t0;
            obj.tf = tf;
            obj.N  = N;
        end
        
        function v = value(obj, x)
            if isa(x, 'yop.ast_independent_initial')
                v = obj(1).evaluate(0);
                
            elseif isa(x, 'yop.ast_independent_final')
                v = obj(end).evaluate(0);
                
            elseif yop.EQ(rem(x-obj.t0, obj.h), 0)
                n = 1 + round((x-obj.t0)/obj.h);
                v = obj(n).evaluate(0);
                
            else
                n = 1 + floor((x-obj.t0)/obj.h);
                tau = rem(x-obj.t0, obj.h)/obj.h;
                v = obj(n).evaluate(tau);
                
            end
        end
        
        function y = vec(obj)
            y = [];
            %for ok = obj(1:end-obj(1).exclude_last)
            for ok = obj(1:end)
                y = [y; ok.y(:)];
            end
        end
        
        function y = mat(obj)
            y = [];
            %for ok = obj(1:end-obj(1).exclude_last)
            for ok = obj(1:end)
                y = [y, ok.y];
            end
        end
        
        function dt = h(obj)
            dt = (obj.tf - obj.t0)/obj.N;
        end
        
    end
end