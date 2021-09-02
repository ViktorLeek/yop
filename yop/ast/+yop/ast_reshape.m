classdef ast_reshape < yop.ast_expression
    properties
        expr
        szs
    end
    methods
        function obj = ast_reshape(expr, varargin)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.szs = varargin;
            obj.dim = size(reshape(ones(size(expr)), varargin{:}));
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            if all(isa_numeric(obj.expr))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function boolv = is_transcription_invariant(obj)
            if all(is_transcription_invariant(obj.expr))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function obj = set_pred(obj)
            add_pred(obj.expr, obj);
        end
        
        function value = evaluate(obj)
            tmp = cell(size(obj.szs));
            for k=1:length(tmp)
                tmp{k} = evaluate(obj.szs{k});
            end
            value = reshape(evaluate(obj.expr), tmp{:});
        end
        
        function v = forward(obj)
            tmp = cell(size(obj.szs));
            for k=1:length(tmp)
                tmp{k} = value(obj.szs{k});
            end
            obj.m_value = reshape(value(obj.expr), tmp{:});
            v = obj.m_value;
        end
        
        function draw(obj)
            str = [];
            for k=1:length(obj.szs)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['reshape(expr, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            draw(obj.expr);
            end_child(obj);
            for k=1:(length(obj.szs)-1)
                begin_child(obj);
                draw(obj.szs{k});
                end_child(obj);
            end
            last_child(obj);
            draw(obj.szs{end});
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
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            for k=1:length(obj.args)
                % probably unnecessary, as args are expected to be numerics
                % but could change in the future.
                [topsort, n_elem, visited] = topological_sort( ...
                    obj.szs{k}, ...
                    visited, ...
                    topsort, ...
                    n_elem ...
                    );
            end
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end