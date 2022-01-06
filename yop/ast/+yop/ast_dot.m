classdef ast_dot < yop.ast_expression
    properties
        A
        B
        d
        nargs
    end
    methods
        function obj = ast_dot(A, B, d)
            obj@yop.ast_expression(is_ival(A) || is_ival(B));
            obj.nargs = nargin;
            switch nargin
                case 2
                    obj.A = A;
                    obj.B = B;
                    obj.dim = size(dot(ones(size(A)), ones(size(B))));
                    
                case 3
                    obj.A = A;
                    obj.B = B;
                    obj.d = d;
                    obj.dim = size(dot(ones(size(A)), ones(size(B)), d));
                    
                otherwise
                    error('Wrong number of arguments provided')
            end
        end
        
        function boolv = isa_numeric(obj)
            boolv = isa_numeric(obj.expr);
            switch obj.nargs
                case 2
                    if all(isa_numeric(obj.A)) &&  all(isa_numeric(obj.B))
                        boolv = true(size(obj));
                    else
                        boolv = false(size(obj));
                    end
                case 3
                    if all(isa_numeric(obj.A)) && all(isa_numeric(obj.B)) && all(isa_numeric(obj.d))
                        boolv = true(size(obj));
                    else
                        boolv = false(size(obj));
                    end
            end
        end
        
        function boolv = isa_reducible(obj)
            if all(isa_reducible(obj.A)) && ...
                    all(isa_reducible(obj.B))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 2
                    value = dot(evaluate(obj.A), evaluate(obj.B));
                    
                case 3
                    value = dot(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d) ...
                        );
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 2
                    obj.m_value = dot(value(obj.A), value(obj.B));
                    
                case 3
                    obj.m_value = dot(...
                        value(obj.A), ...
                        value(obj.B), ...
                        value(obj.d) ...
                        );
            end            
            v = obj.m_value;
        end
        
        
        function draw(obj)
            switch obj.nargs
                case 2
                    fprintf('dot(A, B)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('dot(A, B, dim)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.d);
                    end_child(obj);
                    
            end
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
                topological_sort(obj.A, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.B, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.d, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end