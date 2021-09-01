classdef node < handle
    % Base class for all yop nodes. Primarily used for the AST class.
    
    properties 
        id
        m_value
        pred = {} % predecessors
        dom = {} % dominators
        %idom = {} % immediate dominators
    end
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        % named stream to avoid clash with method.
        stream = yop.stream_state() 
    end
    
    methods
        
        function obj = node()
            obj.id = yop.node.get_uid();
        end
        
        function obj = add_pred(obj, node)
            if isa(obj, 'yop.node')
                obj.pred{end+1} = node;
            end
        end
        
        function obj = reset_pred(obj)
            obj.pred = {};
        end
        
        function i = get_id(obj)
            i = obj.id;
        end
        
        function obj = comp_dom(obj)
            % COMP_DOM - Compute dominators
            
            if isempty(obj.pred)
                obj.dom = {obj};
                return;
            end
            
            % Definition: dom(n0) = {n0}
            %             dom(n) = {n} U {ISEC_{p in pred(n)} dom(p)}
            ds = obj.pred{1}.dom;
            ids = yop.get_ids(ds);
            for k=2:length(obj.pred)
                dk = obj.pred{k}.dom;
                idk = yop.get_ids(dk);
                [ids, idx, ~] = intersect(ids, idk);
                ds = ds(idx);
            end
            obj.dom = {obj, ds{:}};
        end
        
%         function s = sdom(obj)
%             if length(obj.dom) == 1
%                 s = {};
%             else
%                 s = {obj.dom{2:end}};
%             end
%         end
        
        function vars = get_variables(obj)
            [tsort, n_elem] = topological_sort(obj);
            vars = {};
            for k=1:n_elem
                if isa(tsort{k}, 'yop.ast_variable')
                    vars = {vars{:}, tsort{k}};
                end
            end
        end
        
        function reset_stream(obj)
            reset(obj.stream);
        end
        
        function indent(obj)
            % Indents the print by implementing the behvaiour:
            % for k=1:obj.stream.indent_level
            %     if obj.stream.branches(k)
            %         fprintf('|');
            %     else
            %         fprintf(' ');
            %     end
            % end
            
            str = repmat(' ', 1, obj.stream.indent_level-1);
            str(obj.stream.branches) = '|';
            fprintf(str);
        end
        
        function indent_more(obj)
            obj.stream.indent_level = obj.stream.indent_level + 2;
        end
        
        function indent_less(obj)
            obj.stream.indent_level = obj.stream.indent_level - 2;
        end
        
        function begin_child(obj)
            indent(obj);
            fprintf('+-');
            obj.stream.branches(obj.stream.indent_level) = true;
            indent_more(obj);
        end
        
        function end_child(obj)
            indent_less(obj)
            obj.stream.branches(obj.stream.indent_level) = false;
        end
        
        function last_child(obj)
            indent(obj);
            fprintf('+-');
            obj.stream.branches(obj.stream.indent_level) = false;
            indent_more(obj);
        end
        
        function v = value(obj)
            v = obj.m_value;
        end
        
    end
    
    methods (Static)
        function id = get_uid()
            persistent ID
            if isempty(ID)
                ID = 1;
            else
                ID = ID + 1;
            end
            id = ID;
        end
    end
    
end