classdef ast_subsasgn < yop.ast_expression
    
    properties
        node
        s
        b
    end
    
    methods
        function obj = ast_subsasgn(node, s, b)
            obj@yop.ast_expression();
            obj.node = node;
            obj.s = s;
            obj.b = b;
            obj.dim = size(node);
        end
        
        function value = evaluate(obj)
            % Subsref are only created if indices are 's' are numerics, and
            % so they can be passed as they are.
            value = subsasgn(evaluate(obj.node), obj.s, evaluate(obj.b));
        end
        
        function v = forward(obj)
            obj.m_value = subsasgn(value(obj.node), obj.s, value(obj.b));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('subsasgn(node, s, b)\n');
            
            begin_child(obj);
            draw(obj.node);
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
            draw(obj.node);
            end_child(obj);
        end
        
        function [topsort, visited, n_elem] = ...
                topological_sort(obj, topsort, visited, n_elem)
            % Topological sort of expression graph by a dfs.
            
            if nargin == 1
                % Start new sort
                visited = [];
                topsort = cell( ...
                    yop.constants().topsort_preallocation_size, 1);
                n_elem = 0;
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited, n_elem] = ...
                topological_sort(obj.node, topsort, visited, n_elem);
            
            [topsort, visited, n_elem] = ...
                topological_sort(obj.s, topsort, visited, n_elem);
            
            [topsort, visited, n_elem] = ...
                topological_sort(obj.b, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end