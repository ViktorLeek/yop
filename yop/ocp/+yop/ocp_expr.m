classdef ocp_expr < handle
    properties        
        type
        ast
        mx
        sym
        fn
    end
    methods
        function obj = ocp_expr(expr, type)
            obj.ast = expr;
            if nargin == 2
                obj.type = type;
            end
            sz = size(expr);
            obj.mx = yop.cx('expr', sz(1), sz(2));
            obj.sym = sym('expr', sz);
        end
        
        function bool = is_ival(obj)
            bool = is_ival(obj.ast);
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
        
        function obj = set_sym(obj)
            for k=1:length(obj)
                obj(k).ast.m_value = obj(k).sym;
            end
        end

        function obj = set_mx(obj)
            for k=1:length(obj)
                obj(k).ast.m_value = obj(k).mx;
            end
        end
        
        function vec = mx_vec(obj)
            vec = [];
            for k=1:length(obj)
                vec = [vec(:); obj(k).mx(:)];
            end
        end
        
        function vec = sym_vec(obj)
            vec = [];
            for k=1:length(obj)
                vec = [vec(:); obj(k).sym(:)];
            end
        end
        
        function bool = is_transcription_invariant(obj)
            bool = all(is_transcription_invariant(obj.ast));
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