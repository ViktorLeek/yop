classdef ast_node < handle
    properties 
        m_id
        m_value
    end
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        % named stream to avoid clash with method.
        m_stream = yop.stream_state() 
    end
    
    methods
        
        function obj = ast_node(value)
            obj.m_id = yop.ast_node.get_uid();
            obj.m_value = value;
        end
        
        function v = value(obj)
            v = obj.m_value;
        end
        
        function id = ID(obj)
            id = obj.m_id;
        end
        
        function reset_stream(obj)
            reset(obj.m_stream);
        end
        
        function indent(obj)
            % Indents the print by implementing the behvaiour:
            % for k=1:obj.m_stream.indent_level
            %     if obj.m_stream.branches(k)
            %         fprintf('|');
            %     else
            %         fprintf(' ');
            %     end
            % end
            
            str = repmat(' ', 1, obj.m_stream.indent_level-1);
            str(obj.m_stream.branches) = '|';
            fprintf(str);
        end
        
        function indent_more(obj)
            obj.m_stream.indent_level = obj.m_stream.indent_level + 2;
        end
        
        function indent_less(obj)
            obj.m_stream.indent_level = obj.m_stream.indent_level - 2;
        end
        
        function begin_child(obj)
            indent(obj);
            fprintf('+-');
            obj.m_stream.branches(obj.m_stream.indent_level) = true;
            indent_more(obj);
        end
        
        function end_child(obj)
            indent_less(obj)
            obj.m_stream.branches(obj.m_stream.indent_level) = false;
        end
        
        function last_child(obj)
            indent(obj);
            fprintf('+-');
            obj.m_stream.branches(obj.m_stream.indent_level) = false;
            indent_more(obj);
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