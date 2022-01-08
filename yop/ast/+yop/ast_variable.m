classdef ast_variable < yop.ast_expression
    
    properties
        m_name
        m_weight = 1 % Scaling: x_s = (x + os)/w
        m_offset = 0 
    end
    
    methods
        
        function obj = ast_variable(name, weight, offset, isreducible, type)
            obj@yop.ast_expression( ...
                yop.cx(name)         , ... value
                nan                  , ... numval
                yop.initial_timepoint, ... t0
                yop.final_timepoint  , ... tf
                false                , ... isder
                isreducible          , ... isreducible
                type                 , ... type
                0                     ... typeid
                );
            obj.m_typeid = obj.m_id; % 'bit of an ugly fix
            obj.m_name = name;
            obj.m_weight = weight;
            obj.m_offset = offset;
        end
        
        function ast(obj)
            fprintf(['[', num2str(obj.m_id), ']:', obj.m_name, '\n']);
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
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end