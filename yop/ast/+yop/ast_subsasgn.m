classdef ast_subsasgn < yop.ast_expression
    
    properties
        m_expr
        m_s
        m_b
    end
    
    methods
        function obj = ast_subsasgn(expr, s, b)
            
            % Casadi is a bit picky, so logical indices must be converted
            % to indices.
            for k=1:length(s.subs)
                if islogical(s.subs{k})
                    s.subs{k} = find(s.subs{k});
                end
            end
            
            [type_e, tid_e] = Type(expr);
            [type_b, tid_b] = Type(b);
            obj@yop.ast_expression( ...
                subsasgn(value(expr)        , s, value(b))        , ...
                subsasgn(numval(expr)       , s, numval(b))       , ...
                subsasgn(get_t0(expr)       , s, get_t0(b))       , ...
                subsasgn(get_tf(expr)       , s, get_tf(b))       , ...
                subsasgn(get_der(expr)      , s, get_det(b))      , ...
                subsasgn(isa_reducible(expr), s, isa_reducible(b)), ...
                subsasgn(type_e             , s, type_b)          , ...
                subsasgn(tid_e              , s, tid_b)            ...
                );
            obj.m_expr = node;
            obj.m_s = s;
            obj.m_b = b;
        end
        
        function ast(obj)
            fprintf('subsasgn(node, s, b)\n');
            
            begin_child(obj);
            ast(obj.m_expr);
            end_child(obj);
            
            % Subs are numeric
            str = [];
            for k=1:length(obj.m_s.subs)
                str = [str, '[', num2str(obj.m_s.subs{k}), '], '];
            end
            begin_child(obj);
            fprintf(['{', str(1:end-2), '}\n']);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_b);
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
                topological_sort(obj.m_expr, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_s, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_b, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end