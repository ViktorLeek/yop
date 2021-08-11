classdef ast_subsasgn < yop.ast_node
    
    properties
        node
        s
        b
    end
    
    methods
        function obj = ast_subsasgn(node, s, b)
            obj.node = node;
            obj.s = s;
            obj.b = b;
            obj.dim = size(node);
        end
    end
    
    methods % Printing
        
        function ast(obj)
            fprintf('subsasgn(node, s, b)\n');
            
            begin_child(obj);
            ast(obj.node);
            end_child(obj);
            
            % Subs are numeric
            str = [];
            for k=1:length(obj.s.subs)
                str = [str, '[', num2str(obj.s.subs{k}), '], '];
            end
            begin_child(obj);
            fprintf(['{', str(1:end-2), '}\n']);
            end_child(obj);
            
            last_child(obj);
            if isa(obj.b, 'yop.node')
                ast(obj.node);
            else
                % Numeric
                fprintf('[');
                fprintf(num2str(obj.b));
                fprintf(']\n');
            end
            end_child(obj);
        end
        
    end
end