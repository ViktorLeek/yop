classdef sim_poly < yop.lagrange_polynomial
        
    properties
        t0
        tf
    end
    methods
        
        function obj = sim_poly(x, y, t0, tf)
            obj@yop.lagrange_polynomial(x, y);
            obj.t0 = t0;
            obj.tf = tf;
        end
        
        function v = value(obj, x)
            if isa(x, 'yop.ast_independent_initial')
                v = obj(1).evaluate(0);
                
            elseif isa(x, 'yop.ast_independent_final')
                v = obj(end).evaluate(0);
                
            else    
                for n=1:length(obj)
                    if obj(n).t0 <= x && x <= obj(n).tf
                        tau = (x - obj(n).t0)/(obj(n).tf-obj(n).t0);
                        v = obj(n).evaluate(tau);
                        return;
                    end
                end
                error(yop.msg.unexpected_error);
            end
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
        
    end
end