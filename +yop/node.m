classdef node < handle & OOL
    % Implements the basic printing behaviour
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        stream = yop.stream_state() % named stream to avoid clash with method.
    end
    
    methods
        
        function obj = node()
        end
    end
    
    methods % AST      
        
    end
    
    methods % printing
        
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
    end
end