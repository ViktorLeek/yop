classdef ocp_expr < handle
    properties        
        ast
        mx
        sym
        fn
        is_hard
    end
    methods
        function obj = ocp_expr(expr, is_hard)
            obj.ast = expr;
            if nargin == 2
                obj.is_hard = is_hard;
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
        
    end
end