classdef ast_subsref < yop.ast_expression
    
    properties
        node
        s
    end
    
    methods
        function obj = ast_subsref(node, s)
            obj@yop.ast_expression();
            obj.node = node;
            obj.s = s;
            obj.dim = size( subsref( ones(size(node)), s ) );
        end
        
        function value = evaluate(obj)
            value = subsref(evaluate(obj.node), obj.s);
        end
        
        function bool = isa_variable(obj)
            bool = isa_variable(obj.node);
        end
        
        function bool = is_differential(obj)
            bool = is_differential(obj.node);
        end
        
%         function bool = isnumeric(obj)
%             % It would be preferable to inspect the subindices and see if
%             % any of those in obj.node isnumeric. That is however on the
%             % wish list for now.
%         end
        
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