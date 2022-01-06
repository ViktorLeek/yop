classdef ast_subsref < yop.ast_expression
    
    properties
        node
        s
    end
    
    methods
        function obj = ast_subsref(node, s)
            
            % Casadi is a bit picky, so logical indices must be converted
            % to indices.
            for k=1:length(s.subs)
                if islogical(s.subs{k})
                    s.subs{k} = find(s.subs{k});
                end
            end
            obj@yop.ast_expression(is_ival(node));
            obj.node = node;
            obj.s = s;
            obj.dim = size( subsref( ones(size(node)), s ) );
        end
        
        function val = numval(obj)
            val = subsref(numval(obj.node), obj.s);           
        end
        
        function [type, id] = Type(obj)
            if ~isempty(obj.m_type)
                type = obj.m_type.type;
                id = obj.m_type.id;
                return;
            end
            [type, id] = Type(obj.node);
            idx = get_indices(obj);
            type = type(idx);
            id = id(idx);
            obj.m_type.type = type;
            obj.m_type.id = id;
        end
        
        function [bool, id] = isa_der(obj)
            [bool, id] = isa_der(obj.node);
            idx = get_indices(obj);
            bool = bool(idx);
            id = id(idx);
        end
        
        function [bool, tp] = isa_timepoint(obj)
            [bool, tp] = isa_timepoint(obj.node);
            idx = get_indices(obj);
            bool = bool(idx);
            tp = tp(idx);
        end
        
        function boolv = isa_reducible(obj)
            boolv = isa_reducible(obj.node);
            boolv = boolv(get_indices(obj));
        end
        
        function idx = get_indices(obj)
            % Return the indices this subsref node refers to
            % persistant variable for idx?
            node_val = obj.node.m_value;
            sz = size(obj.node);
            obj.node.m_value = reshape(1:prod(sz), sz);
            idx = subsref(value(obj.node), obj.s);
            idx = idx(:);
            obj.node.m_value = node_val;
        end
        
        function value = evaluate(obj)
            value = subsref(evaluate(obj.node), obj.s);
        end
        
        function v = forward(obj)
            obj.m_value = subsref(value(obj.node), obj.s);
            v = obj.m_value;
        end
        
        function ast(obj)
            fprintf('subsref(node, s)\n');
            
            begin_child(obj);
            ast(obj.node);
            end_child(obj);
            
            str = [];
            for k=1:length(obj.s.subs)
                idx = obj.s.subs{k};
                if size(idx,1) > size(idx,2)
                    idx = idx';
                end
                str = [str, '[', num2str(idx), '], '];
            end
            last_child(obj);
            fprintf(['{', str(1:end-2), '}\n']);
            end_child(obj);
        end
        
        function [topsort, n_elem, visited] = ...
                topological_sort(obj, visited, topsort, n_elem)
            % Topological sort of expression graph by a dfs.
            
            switch nargin
                case 1
                    % Start new sort: topological_sort(obj)
                    visited = [];
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size, 1);
                    n_elem = 0;
                    
                case 2
                    % Semi-warm start: topological_sort(obj, visited)
                    % In semi-warm start 'visited' is already provided, but 
                    % no elements are sorted. This is for instance useful 
                    % for finding all variables in a number of expressions 
                    % that are suspected to contain common subexpressions.
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size, 1);
                    n_elem = 0;
                    
                otherwise
                    % Pass
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.node, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.s, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end