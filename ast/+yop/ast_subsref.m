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
        
        function v = forward(obj)
            obj.m_value = subsref(value(obj.node), obj.s);
            v = obj.m_value;
        end
        
        function bool = isa_variable(obj)
            bool = isa_variable(obj.node);
            
            % The bools vector is for every element in obj.node. Since this
            % node only refers to a few of those elements, it is necessary
            % to extract those. This is done by enumerating all of the
            % elements of obj.node, forward evaluating this node, and then
            % the indices are are the value of this node. 
            sz = size(obj.node);
            obj.node.m_value = reshape(1:prod(sz), sz);
            idx_matrix = forward(obj);
            bool = bool(idx_matrix(:)); 
            
        end
        
        function var = get_variable(obj)
            var  = get_variable(obj.node);
        end
        
        function bool = is_differential(obj)
            bool = is_differential(obj.node);
        end
        
        function bool = is_algebraic(obj)
            bool = is_algebraic(obj.node);
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
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;
            
        end
        
    end
end