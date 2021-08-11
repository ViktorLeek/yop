classdef ast_subsref < yop.ast_node
    
    properties
        node
        s
    end
    
    methods
        function obj = ast_subsref(node, s)
            obj.node = node;
            obj.s = s;
            obj.dim = size( subsref( ones(size(node)), s ) );
        end
    end
    
    methods % Printing the ast
        
        function ast(obj)
            fprintf('subsref(node, s)\n');
            
            begin_child(obj);
            ast(obj.node);
            end_child(obj);
            
            
            str = [];
            for k=1:length(obj.s.subs)
                str = [str, '[', num2str(obj.s.subs{k}), '], '];
            end
            last_child(obj);
            fprintf(['{', str(1:end-2), '}\n']);
            end_child(obj);
        end
        
    end
end