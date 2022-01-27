classdef ast_logical < yop.ast_node
    
    properties
        m_lhs
        m_rhs
    end
    
    methods
        function obj = ast_logical(value, lhs, rhs)
            obj@yop.ast_node(value);
            obj.m_lhs = lhs;
            obj.m_rhs = rhs;
        end
        
        function bool = is_hard(obj)
            bool = false(size(obj));
        end
        
        function bool = isa_reducible(obj)
            bool = false(size(obj));
        end
        
        function sz = size(obj, varargin)
            sz = size(obj.m_value, varargin{:});
        end
        
        function rel = or(lhs, rhs)
            rel = yop.ast_or(lhs, rhs);
        end
        
        function rel = and(lhs, rhs)
            rel = yop.ast_and(lhs, rhs);
        end
        
        function rel = not(lhs, rhs)
            rel = yop.ast_not(lhs, rhs);
        end
        
        function node = if_else(varargin)
            node = yop.ast_if_else(varargin{:});
        end
        
        function idx = end(obj, k, n)
            idx = builtin('end', ones(size(obj)), k, n);
        end
        
        function rels = get_relations(obj)
            l = get_relations(obj.m_lhs);
            r = get_relations(obj.m_rhs);
            rels = {obj, l{:}, r{:}};
            rels = rels(~cellfun('isempty', rels));
        end
        
        function r = rmost(obj)
            r = rmost(obj.m_rhs);
        end
        
        function l = lmost(obj)
            l = lmost(obj.m_lhs);
        end
        
        function ast(obj)
            fprintf([obj.m_name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            ast(obj.m_lhs);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_rhs);
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
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_lhs, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_rhs, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;
            
        end
        
    end
end