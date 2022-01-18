classdef interpolating_poly < yop.lagrange_polynomial
    properties
        t0
        tf
    end
    
    methods
        function obj = interpolating_poly(x, y, t0, tf)
            obj@yop.lagrange_polynomial(x, y);
            obj.t0 = t0;
            obj.tf = tf;
        end
        
        function v = eval(obj, tau)
            % lazy evaluation
            v = [];
            
            if length(obj.x) == 1
                v = obj.y(:);
                return;
            end
            
            [M, I] = min(abs(tau - obj.x));
            if yop.EQ(M, 0, 1e-9) && ~isempty(obj.y)
                % if tau is a collocation point, the call to evaluate is
                % avoided as it is slow.
                v = obj.y(:,I);
            else
                v = obj.evaluate(tau);
            end
        end
        
        function v = value(obj, x)
            % Value at timepoint x
            % obj is a vector or polynomials
            if x == yop.initial_timepoint
                v = obj(1).evaluate(0);
                
            elseif x == yop.final_timepoint && length(obj(end).x)==1
                v = obj(end).evaluate(0);
                
            elseif x == yop.final_timepoint
                v = obj(end).evaluate(1); % assuming x in [0, 1]
                
            elseif yop.EQ(rem(x-obj(1).t0, obj.h), 0)
                n = 1 + round((x-obj(1).t0)/obj.h);
                if n > length(obj)
                    v = obj(end).evaluate(1);
                else
                    v = obj(n).evaluate(0);
                end
                
            else
                n = 1 + floor((x-obj(1).t0)/obj.h);
                tau = rem(x-obj(1).t0, obj.h)/obj.h;
                v = obj(n).evaluate(tau);
                
            end
        end
        
        function v = valuenh(obj, x)
            % Value at timepoint x for non-homogenuous grids
            % obj is a vector or polynomials
            if x == yop.initial_timepoint
                v = obj(1).evaluate(0);
                return
            end
                
            if x == yop.final_timepoint && length(obj(end).x)==1
                v = obj(end).evaluate(0);
                return
            end
                
            if x == yop.final_timepoint
                v = obj(end).evaluate(1); % assuming x in [0, 1]
                return
            end
            
            for o=obj
                if o.t0 <= x && x <= o.tf
                    tau = (x-o.t0)/o.dt;
                    v = o.evaluate(tau);
                    return
                end
            end
            
            error(yop.error.timepoint_outside_range());
        end
        
        function y = vec(obj)
            y = [];
            for ok = obj(1:end)
                y = [y; ok.y(:)];
            end
        end
        
        function y = mat(obj)
            y = [];
            for ok = obj(1:end)
                y = [y, ok.y];
            end
        end
        
        function dt = h(obj)
            % Uniform grid
            dt = (obj(end).tf - obj(1).t0)/obj.N;
        end
        
        function h = dt(obj)
            % Non-uniform grid
            h = obj.tf - obj.t0;
        end
        
        function n = N(obj)
            % Number of segments
            n = length(obj);
            if length(obj(end).x) == 1
                n = n - 1;
            end
        end
        
    end
end