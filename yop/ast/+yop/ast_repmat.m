classdef ast_repmat < yop.ast_expression
    properties
        m_expr
        m_args
    end
    methods
        function obj = ast_repmat(expr, varargin)
            obj@yop.ast_expression( ...
                repmat(expr.m_value    , varargin{:}), ... value
                repmat(expr.m_numval   , varargin{:}), ... numval
                repmat(expr.m_t0       , varargin{:}), ... t0
                repmat(expr.m_tf       , varargin{:}), ... tf
                repmat(expr.m_der      , varargin{:}), ... der
                repmat(expr.m_reducible, varargin{:}), ... reducible
                repmat(expr.m_type     , varargin{:}), ... type
                repmat(expr.m_typeid   , varargin{:}) ... typeid
                );
            obj.m_expr = expr;
            obj.m_args = varargin;
        end
        
        function ast(obj)
            str = [];
            for k=1:length(obj.m_args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['repmat(expr, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.m_expr);
            end_child(obj);
            for k=1:(length(obj.m_args)-1)
                begin_child(obj);
                ast(obj.m_args{k});
                end_child(obj);
            end
            last_child(obj);
            ast(obj.m_args{end});
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
            if ~isempty( find(visited == obj.m_id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.m_id];
            
            % Visit child
            for k=1:length(obj.m_args)
                % probably unnecessary, as args are expected to be numerics
                % but could change in the future.
                [topsort, n_elem, visited] = topological_sort( ...
                    obj.m_args{k}, ...
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