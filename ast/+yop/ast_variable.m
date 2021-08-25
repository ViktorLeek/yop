classdef ast_variable < yop.ast_expression
    % The meaning of an ast_variable is a problem variable in the OCP.
    % Variables and expressions have almost the same semantics, but they
    % differ in a few ways:
    %   1) Variables can represent time-varying and constant values,
    %      whereas expressions are always considered time-varying. This
    %      means that states for instance can be evaluated at timepoints,
    %      whereas parameters cannot.
    %   2) Time-varying variables can only be eval
    
    properties
        name
    end
    
    methods
        
        function obj = ast_variable(name, rows, cols)    
            obj@yop.ast_expression();
            obj.name = name;
            obj.dim = [rows, cols];
        end
        
        function bool = isa_variable(obj)
            bool = true(size(obj));
        end
        
        function boolv = isa_numeric(obj)
            boolv = false(size(obj));
        end
        
        function draw(obj)
            fprintf([obj.name, '\n']);
        end
        
        function value = evaluate(obj)
            value = obj.value;
        end
        
        function v = forward(obj)
            v = obj.m_value;
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
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end