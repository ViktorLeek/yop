classdef ocp_int < handle & matlab.mixin.Copyable
    properties        
        node
        mx
        sym
        fn
        disc
    end
    methods
        function obj = ocp_int(ast_int)
            obj.node = ast_int;
        end
        
        function expr = expr(obj)
            expr = obj.node.expr;
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
        
%         function vec = get_disc(obj)
%             vec = [];
%             for k=1:length(obj)
%                 if isempty(obj(k).disc)
%                     tmp = zeros(size(obj(k).expr));
%                     vec = [vec(:); tmp(:)];
%                 else
%                     vec = [vec(:); obj(k).disc];
%                 end
%             end
%         end
        
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