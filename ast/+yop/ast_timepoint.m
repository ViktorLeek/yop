classdef ast_timepoint < yop.ast_expression
    properties
        timepoint
        expr
    end
    methods
        function obj = ast_timepoint(tp, expr)
            obj@yop.ast_expression();
            
            if isa(expr, 'yop.ast_timepoint')
                error(['[yop] Error: Cannot evaluate at timepoint at a' ...
                    ' timepoint']);
            end
            
            switch class(tp)
                case 'yop.ast_eq'
                    if isa(tp.lhs, 'yop.ast_independent') && ...
                            isnumeric(tp.rhs)
                        % t == 2
                        timepoint = tp.rhs;
                        
                    elseif isnumeric(tp.lhs) && ...
                            isa(tp.rhs, 'yop.ast_independent')
                        % 2 == t
                        timepoint = tp.lhs;
                        
                    elseif isa(tp.lhs, 'yop.ast_independent') && ...
                            (isa(tp.rhs, 'yop.ast_independent_initial')...
                            || isa(tp.rhs, 'yop.ast_independent_final'))
                         % t == t0 || t == tf
                         timepoint = tp.rhs;
                         
                    elseif (isa(tp.lhs, 'yop.ast_independent_initial')...
                           || isa(tp.lhs, 'yop.ast_independent_final'))...
                            && isa(tp.rhs, 'yop.ast_independent')
                        % t0 == t || tf == t
                        timepoint = tp.lhs;
                        
                    else
                        m='[yop] Error: Illegal relation for a timepoint';
                        error(m);
                        
                    end
                    
                case {'yop.ast_independent', ...
                      'yop.ast_independent_initial' ...
                      'yop.ast_independent_final'}
                  timepoint = tp;
                  
                otherwise
                    error('[yop] Error: Illegal relation for a timepoint');
            end
            
            obj.timepoint = timepoint;
            obj.expr = expr;
        end
        
        function draw(obj)
            fprintf('timepoint(timepoint, expr)\n');
            
            begin_child(obj);
            draw(obj.timepoint);
            end_child(obj);
            
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
        
        function [topsort, visited, n_elem] = ...
                topological_sort(obj, topsort, visited, n_elem)
            % Topological sort of expression graph by a dfs.
            
            if nargin == 1
                % Start new sort
                visited = [];
                topsort = cell( ...
                    yop.constants().topsort_preallocation_size, 1);
                n_elem = 0;
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited, n_elem] = ...
                topological_sort(obj.timepoint, topsort, visited, n_elem);
            
            [topsort, visited, n_elem] = ...
                topological_sort(obj.expr, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end