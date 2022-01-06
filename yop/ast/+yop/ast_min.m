classdef ast_min < yop.ast_expression
    properties
        A
        B
        d
        flag
        nargs
    end
    methods
        function obj = ast_min(A, B, d, flag)
            obj@yop.ast_expression(is_ival(A) || is_ival(B));
            obj.A = A;
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.dim = size(min(ones(size(A))));
                    
                case 2
                    obj.B = B;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp));
                    
                case 3
                    obj.B = B;
                    obj.d = d;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d));
                    
                case 4
                    obj.B = B;
                    obj.d = d;
                    obj.flag = flag;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d, flag));     
            end
            
        end
        
        function val = numval(obj)
            switch obj.nargs
                case 1
                    val = min(numval(obj.A));
                    
                case 2
                    val = min(numval(obj.A), numval(obj.B));
                    
                case 3
                    val = min(...
                        numval(obj.A), ...
                        numval(obj.B), ...
                        numval(obj.d) ...
                        );
                    
                case 4
                    val = min(...
                        numval(obj.A), ...
                        numval(obj.B), ...
                        numval(obj.d), ...
                        numval(obj.flag) ...
                        );
            end
        end
        
        function boolv = isa_reducible(obj)
            switch obj.nargs
                case 1
                    if all(isa_reducible(obj.A))
                        boolv = true(size(obj));
                    else
                        boolv = false(size(obj));
                    end
                    
                case {2, 3, 4}
                    if all(isa_reducible(obj.A)) && ...
                            all(isa_reducible(obj.B))
                        boolv = true(size(obj));
                    else
                        boolv = false(size(obj));
                    end
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = min(evaluate(obj.A));
                    
                case 2
                    value = min(evaluate(obj.A), evaluate(obj.B));
                    
                case 3
                    value = min(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d) ...
                        );
                    
                case 4
                    value = min(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d), ...
                        evaluate(obj.flag) ...
                        );
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 1
                    obj.m_value = min(value(obj.A));
                    
                case 2
                    obj.m_value = min(value(obj.A), value(obj.B));
                    
                case 3
                    obj.m_value = min(...
                        value(obj.A), ...
                        value(obj.B), ...
                        value(obj.d) ...
                        );
                    
                case 4
                    obj.m_value = min(...
                        value(obj.A), ...
                        value(obj.B), ...
                        value(obj.d), ...
                        value(obj.flag) ...
                        );
            end         
            v = obj.m_value;
        end
        
        function draw(obj)
            switch obj.nargs
                case 1
                    fprintf('min(A)\n');
                    last_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('min(A, B)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('min(A, B, dim)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.d);
                    end_child(obj);
                    
                case 4
                    fprintf('min(A, B, dim, flag)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.d);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.flag);
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
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.flag, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end