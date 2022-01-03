classdef interpolating_poly < yop.lagrange_polynomial
    properties
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
            v = [];
            if x == yop.initial_timepoint
                v = obj(1).evaluate(0);
                
            elseif x == yop.final_timepoint && length(obj(end).x)==1
                v = obj(end).evaluate(0);
                
            elseif x == yop.final_timepoint
                v = obj(end).evaluate(1); % assuming x in [0, 1]
                
            elseif yop.EQ(rem(x-obj(1).t0, obj.h), 0)
                n = 1 + round((x-obj(1).t0)/obj.h);
                if n > length(obj)
                    obj(end).evaluate(1);
                else
                    v = obj(n).evaluate(0);
                end
                
            else
                n = 1 + floor((x-obj(1).t0)/obj.h);
                tau = rem(x-obj(1).t0, obj.h)/obj.h;
                v = obj(n).evaluate(tau);
                
            end
        end
        
        %         function [n, r] = get_next_index(obj, x)
        %             % Compute the index to the first value larger than x.
        %             if yop.EQ(rem(x-obj(1).t0, obj.h), 0)
        %                 n = 1 + round((x-obj(1).t0)/obj.h);
        %                 n = max(n, length(obj));
        %                 r = 1;
        %             else
        %                 n = 1 + floor((x-obj(1).t0)/obj.h);
        %                 tau = rem(x-obj(1).t0, obj.h)/obj.h;
        %                 tau_vec = obj(1).x;
        %                 dtau = tau_vec - tau;
        %                 r = find(dtau >= 0, 1);
        %                 if isempty(r)
        %                     % Closest point is next interval.
        %                     n = max(n+1, length(obj));
        %                     r = 1;
        %                 end
        %             end
        %         end
        %
        %         function [n, r] = get_prev_index(obj, x)
        %             % Compute the index to the first value larger than x.
        %             if yop.EQ(rem(x-obj(1).t0, obj.h), 0)
        %                 n = 1 + round((x-obj(1).t0)/obj.h);
        %                 n = max(n, length(obj));
        %                 r = 1;
        %             else
        %                 n = 1 + floor((x-obj(1).t0)/obj.h);
        %                 tau = rem(x-obj(1).t0, obj.h)/obj.h;
        %                 tau_vec = obj(1).x;
        %                 dtau = tau_vec - tau;
        %                 r = find(dtau <= 0, 1, 'last');
        %                 if isempty(r)
        %                     r = 1;
        %                 end
        %             end
        %         end
        
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
            dt = (obj(1).tf - obj(1).t0)/obj(1).N;
        end
        
    end
end