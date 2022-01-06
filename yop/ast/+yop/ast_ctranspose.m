classdef ast_ctranspose < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_ctranspose(expr)
            obj@yop.ast_expression(is_ival(expr));
            obj.expr = expr;
            obj.dim = size(ctranspose(ones(size(expr))));
        end
        
        function val = numval(obj)
            val = ctranspose(numval(obj.expr));
        end
        
        function [type, id] = Type(obj)
            if ~isempty(obj.m_type)
                type = obj.m_type.type;
                id = obj.m_type.id;
                return;
            end
            [type, id] = Type(obj.expr);
            type = ctranspose(type);
            id = ctranspose(id);
            obj.m_type.type = type;
            obj.m_type.id = id;
        end
        
        function [bool, tp, type] = isa_timepoint(obj)
            [bool, tp, type] = isa_timepoint(obj.expr);
            bool = ctranspose(bool);
            tp = ctranspose(tp);
            type = ctranspose(type);
        end
        
        function [bool, id] = isa_der(obj)
            [bool, id] = isa_der(obj.expr);
            bool = ctranspose(bool);
            id = ctranspose(id);
        end
        
        function boolv = isa_reducible(obj)
            boolv = isa_reducible(obj.expr);
        end
        
        function value = evaluate(obj)
            value = ctranspose(evaluate(obj.expr));
        end
        
        function v = forward(obj)
            obj.m_value = ctranspose(value(obj.expr));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('ctranspose(obj)\n');
            last_child(obj);
            draw(obj.expr);
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
                        cell(yop.constants().topsort_preallocation_size,1);
                    n_elem = 0;
                    
                case 2
                    % Semi-warm start: topological_sort(obj, visited)
                    % In semi-warm start 'visited' is already provided, but 
                    % no elements are sorted. This is for instance useful 
                    % for finding all variables in a number of expressions 
                    % that are suspected to contain common subexpressions.
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size,1);
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
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end