classdef node < handle
    % Base class for all yop nodes. Primarily used for the AST class.
    
    properties 
        id
        m_value
    end
    
    %     properties %(to remove?)
    %         pred = {} % predecessors
    %         dom = {} % dominators
    %         %idom = {} % immediate dominators
    %     end
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        % named stream to avoid clash with method.
        stream = yop.stream_state() 
    end
    
    methods
        
        function obj = node()
            obj.id = yop.node.get_uid();
        end
        
        function node = at(expression, timepoint)
            % Alternative syntax for evaluating expression at a timepoint
            node = yop.ast_timepoint(timepoint, expression);
        end
        
        function value = fw_eval(expr)
            % FW_EVAL - Forward evaluate
            [sort, K] = topological_sort(expr);
            for k=1:(K-1)
                forward(sort{k});
            end
            value = forward(sort{K});
        end
        
        function val = propagate_value(expr)
            [sort, K] = topological_sort(expr);
            for k=1:K
                switch class(sort{k})
                    case {'yop.ast_int', 'yop.ast_timepoint', 'yop.ast_der'}
                        % These do not propagate values, so we need a 
                        % manual bridge here.
                        sort{k}.m_value = value(sort{k}.expr);
                    otherwise
                        forward(sort{k});
                end
            end
            val = value(sort{K});
        end
        
        
        
        function id = get_id(obj)
            id = obj.id;
        end
        
        function vars = get_vars(obj)
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