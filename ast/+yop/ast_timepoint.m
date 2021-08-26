classdef ast_timepoint < yop.ast_expression
    properties
        timepoint
        expr
    end
    methods
        function obj = ast_timepoint(tp, expr)
            obj@yop.ast_expression();
            obj.dim = expr.dim;
            if isa(expr, 'yop.ast_timepoint')
                error(['[yop] Error: Cannot evaluate at timepoint at a' ...
                    ' timepoint']);
            end
            
            switch class(tp)
                case 'yop.ast_eq'
                    if isa(tp.lhs, 'yop.ast_independent') && ...
                            isnumeric(tp.rhs)
                        % t == 2
                        if isinf(tp.rhs)
                            error(['[Yop] Error: inf is not a '...
                                'valid timepoint']);
                        end
                        timepoint = tp.rhs;
                        
                    elseif isnumeric(tp.lhs) && ...
                            isa(tp.rhs, 'yop.ast_independent')
                        % 2 == t
                        if isinf(tp.lhs)
                            error(['[Yop] Error: inf is not a '...
                                'valid timepoint']);
                        end
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
                    
                case {'yop.ast_independent_initial', ...
                      'yop.ast_independent_final'}
                  timepoint = tp;
                  
                case 'yop.ast_independent'
                  error(['[yop] Error: Do not specify the independent', ...
                      ' variable as a timepoint. It represents the '...
                      'evolution of time, not a timepoint'])
                  
                otherwise
                    error('[yop] Error: Illegal relation for a timepoint');
            end
            
            obj.timepoint = timepoint;
            obj.expr = expr;
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            boolv = isa_numeric(obj.expr);
        end
        
        function [bool, tp] = isa_timepoint(obj)
            bool = true(size(obj.expr));
            switch class(obj.timepoint)    
                case 'yop.ast_independent_initial'
                    tp = yop.initial_timepoint(size(obj.expr));
                    
                case 'yop.ast_independent_final'
                    tp = yop.final_timepoint(size(obj.expr));
                    
                otherwise
                    % numeric types
                    tp = obj.timepoint*ones(size(obj.expr));
            end
        end
        
        function value = evaluate(obj)
            % This simply propagates through. The purpose is to be able to
            % go to single variable form by enumaration, so this function
            % must be able to propagate indices.
            value = evaluate(obj.expr);
        end
        
        function v = forward(obj)
            % This simply propagates through. The purpose is to be able to
            % go to single variable form by enumaration, so this function
            % must be able to propagate indices.
            obj.m_value = value(obj.expr);
            v = obj.m_value;
        end
        
        function [bool, id] = isa_variable(obj)
            [bool, id] = isa_variable(obj.expr);
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
                topological_sort(obj.timepoint, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end