classdef ocp_expr < handle
    properties        
        type
        ast
        mx
        sym
        fn
        is_hard
    end
    methods
        function obj = ocp_expr(expr, aux)
            obj.ast = expr;
            if nargin == 2
                if isnumeric(aux)
                    obj.type = aux;
                else
                    obj.is_hard = aux;
                end
            end
            sz = size(expr);
            obj.mx = casadi.MX.sym('expr', sz(1), sz(2));
            obj.sym = sym('expr', sz);
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

%         function [sn, n_tp, n_int, n_der] = mx_functionx(obj, t, x, u, p) 
%             % x for extra parameters
% 
%             [sn, tps, ints, ders] = yop.ocp.find_special_nodes(obj.ast);
%             
%             n_tp = n_elem(tps);
%             n_int = n_elem(ints);
%             n_der = n_elem(ders);
%             
%             set_mx([t, x, u, p]);
%             set_mx([tps, ints, ders]);
%             
%             args =  {mx_vec(t), mx_vec(x), mx_vec(u), mx_vec(p), ...
%                 mx_vec(tps), mx_vec(ints)};
%             
%             for node = [tps, ints]
%                 mx_expr = fw_eval(node.ast.expr);
%                 node.fn = casadi.Function('fn', args, {mx_expr});
%             end
%             
%             mx_expr = fw_eval(obj.ast);
%             obj.fn = casadi.Function('fn', args, {mx_expr});
%         end
%         
%         function mx_function(obj, t, x, u, p)
%             set_mx([t, x, u, p]);
%             mx_expr = fw_eval(obj.ast);
%             args = {mx_vec(t), mx_vec(x), mx_vec(u), mx_vec(p)};
%             obj.fn = casadi.Function('fn', args, {mx_expr});
%         end