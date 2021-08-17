classdef ast_variable < yop.ast_expression
    
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
            bool = true;
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
        
        function [topsort, visited] = topological_sort(obj, topsort, visited)
            % Topological sort of expression graph by a dfs.
            
            % Initialize if second and third args are empty
            if nargin == 1
                topsort = {};
                visited = [];
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
        
    end
end