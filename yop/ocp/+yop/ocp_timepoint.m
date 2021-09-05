classdef ocp_timepoint < handle
    properties        
        node
        mx
        sym
        fn
        disc
    end
    methods
        function obj = ocp_timepoint(ast_timepoint)
            obj.node = ast_timepoint;
        end
        
        function expr = expr(obj)
            expr = obj.node.expr;
        end
        
        function tp = timepoint(obj)
            tp = obj.node.timepoint;
        end
        
        function obj = set_sym(obj)
            for k=1:length(obj)
                obj(k).node.m_value = obj(k).sym;
            end
        end

        function obj = set_mx(obj)
            for k=1:length(obj)
                obj(k).node.m_value = obj(k).mx;
            end
        end
        
        function vec = get_disc(obj)
            vec = [];
            for k=1:length(obj)
                if isempty(obj(k).disc)
                    tmp = zeros(size(obj(k).expr));
                    vec = [vec(:); tmp(:)];
                else
                    vec = [vec(:); obj(k).disc];
                end
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