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
            obj.t0 = t0;
            obj.tf = tf;
            obj.expr = expr;
        end
        
        function val = numval(obj)
            val = numval(obj.expr);
        end
        
        function [type, id] = Type(obj)
            if ~isempty(obj.m_type)
                type = obj.m_type.type;
                id = obj.m_type.id;
                return;
            end
            [type, id] = Type(obj.expr);
            obj.m_type.type = type;
            obj.m_type.id = id;
        end
        
        function boolv = isa_reducible(obj)
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
        
        function ast(obj)
            fprintf(['[', num2str(obj.id), ']:', ...
                'interval(t0, tf, expr)\n']);
            
            begin_child(obj);
            ast(obj.t0);
            end_child(obj);
            
            begin_child(obj);
            ast(obj.tf);
            end_child(obj);
            
            last_child(obj);
            ast(obj.expr);
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