classdef node < handle
    % Base class for all yop nodes. Primarily used for the AST class.
    
    properties 
        id
        m_value
    end
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        % named stream to avoid clash with method.
        stream = yop.stream_state() 
    end
    
    methods
        
        function obj = node()
            obj.id = yop.node.get_id();
        end
        
        function reset_stream(obj)
            reset(obj.stream);
        end
        
        function vars = get_variables(obj)
            [tsort, n_elem] = topological_sort(obj);
            vars = {};
            for k=1:n_elem
                if isa(tsort{k}, 'yop.ast_variable')
                    vars = {vars{:}, tsort{k}};
                end
            end
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
        function id = get_id()
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