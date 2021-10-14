classdef ast_timeinterval < yop.ast_expression
    properties
        t0
        tf
        expr
    end
    methods
        function obj = ast_timeinterval(t0, tf, expr)
            obj@yop.ast_expression(true);
            obj.dim = expr.dim;
            
            if isnumeric(t0) 
                obj.t0 = t0;
            elseif isa(t0, 'yop.ast_independent_initial')
                obj.t0 = yop.initial_timepoint;
            else
                error(yop.msg.illegal_timepoint);
            end
            
            if isnumeric(tf) 
                obj.tf = tf;
            elseif isa(tf, 'yop.ast_independent_final')
                obj.tf = yop.final_timepoint;
            else
                error(yop.msg.illegal_timepoint);
            end
            
            obj.expr = expr;
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            boolv = isa_numeric(obj.expr);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = false(size(obj));
        end
        
        function value = evaluate(obj)            
            % Evaluates like an ast_variable. The reason is that it does
            % not have the semantics of a timevarying expression, so it
            % needs to be parameterized separately. Before parametrization
            % has occured, the idea is that an MX/sym variable (of the
            % right size) is used to represent its value.
            value = evaluate(obj.expr);
        end
        
        function v = forward(obj)
            % Evaluates like an ast_variable. The reason is that it does
            % not have the semantics of a timevarying expression, so it
            % needs to be parameterized separately. Before parametrization
            % has occured, the idea is that an MX/sym variable (of the
            % right size) is used to represent its value.
            obj.m_value = value(obj.expr);
            v = obj.m_value;
        end
        
        function [bool, id] = isa_variable(obj)
            [bool, id] = isa_variable(obj.expr);
        end
        
        function [bool, id] = isa_state(obj)
            [bool, id] = isa_state(obj.expr);
        end
        
        function [bool, id] = isa_independent(obj)
            [bool, id] = isa_independent(obj.expr);
        end
        
        function [bool, id] = isa_parameter(obj)
            [bool, id] = isa_parameter(obj.expr);
        end
        
        function [bool, id] = isa_control(obj)
            [bool, id] = isa_control(obj.expr);
        end
        
        function [bool, id] = isa_algebraic(obj)
            [bool, id] = isa_algebraic(obj.expr);
        end
        
        function draw(obj)
            fprintf(['[', num2str(obj.id), ']:', ...
                'interval(timepoint, expr)\n']);
            
            begin_child(obj);
            draw(obj.ival);
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
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end