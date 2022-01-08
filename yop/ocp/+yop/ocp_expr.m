classdef ocp_expr < handle
    properties        
        type
        ast
        fn
    end
    methods
        function obj = ocp_expr(expr, type)
            obj.ast = expr;
            if nargin == 2
                obj.type = type;
            end
        end
        
        function bool = isa_ival(obj)
            bool = isa_ival(obj.ast);
        end
        
        
        function [t0, tf] = get_ival(obj)
            [t0, tf] = get_ival(obj.ast);
        end
        
        function n = n_elem(obj)
            n = 0;
            for k=1:length(obj)
                n = n + prod(size(obj(k).ast));
            end
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for o=obj
                vec = [vec(:); o.ast.m_value(:)];
            end
        end
        
        function bool = isa_reducible(obj)
            bool = all(isa_reducible(obj.ast));
        end
        
        function tp = timepoint(obj)
            % Notice, error if not a timepoint.
            tp = obj.ast.timepoint;
        end
        
    end
    methods (Static)
        function t = tp()
            t = 1;
        end
        function t = int()
            t = 2;
        end
        function t = der()
            t = 3;
        end
    end
end