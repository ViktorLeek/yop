classdef lagrange_signal < yop.lagrange_polynomial
    properties
        exclude_last = false;
    end
    
    methods
        function obj = lagrange_signal(x, y)
            obj@yop.lagrange_polynomial(x, y);
        end
        
        function v = value(obj, x, t0, tf)
            if isa(x, 'yop.ast_independent_initial')
                v = obj(1).evaluate(0);
            elseif isa(x, 'yop.ast_independent_final')
                v = obj(end).evaluate(0);
            else    
                h = (tf - t0)/(length(obj)-1);
                n = 1 + floor((x-t0)/h);
                tau = rem(x-t0, h)/h;
                v = obj(n).evaluate(tau);
            end
        end
        
        function y = vec(obj)
            y = [];
            for ok = obj(1:end-obj(1).exclude_last)
                y = [y; ok.y(:)];
            end
        end
        
        function y = mat(obj)
            y = [];
            for ok = obj(1:end-obj(1).exclude_last)
                y = [y, ok.y];
            end
        end
        
        %         function t = t0(obj)
        %             t = obj(1).T0;
        %         end
        %
        %         function dt = h(obj)
        %             % last object is a point, therefore -1.
        %             dt = (obj(1).tf - obj(1).t0)/(length(obj)-1);
        %         end
    end
end