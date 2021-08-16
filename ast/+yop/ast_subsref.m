classdef ast_subsref < yop.ast_expression
    
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
        
        function value = evaluate(obj)
            value = subsref(evaluate(obj.node), obj.s);
        end
        
        function draw(obj)
            fprintf('subsref(node, s)\n');
            
            begin_child(obj);
            draw(obj.node);
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