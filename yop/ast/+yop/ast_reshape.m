classdef ast_reshape < yop.ast_expression
    properties
        expr
        szs
    end
    methods
        function obj = ast_reshape(expr, varargin)
            obj@yop.ast_expression(is_ival(expr));
            obj.expr = expr;
            obj.szs = varargin;
            obj.dim = size(reshape(ones(size(expr)), varargin{:}));
        end
        
        function val = numval(obj)
            val = reshape(numval(obj.expr), obj.szs{:});
        end
        
        function [type, id] = Type(obj)
            if ~isempty(obj.m_type)
                type = obj.m_type.type;
                id = obj.m_type.id;
                return;
            end
            [type, id] = Type(obj.expr);
            type = reshape(type, obj.szs{:});
            id = reshape(id, obj.szs{:});
            obj.m_type.type = type;
            obj.m_type.id = id;
        end
        
        function [bool, tp] = isa_timepoint(obj)
            [bool, tp] = isa_timepoint(obj.expr);
            bool = reshape(bool, obj.szs{:});
            tp = reshape(tp, obj.szs{:});
        end
        
        function [bool, id] = isa_der(obj)
            [bool, id] = isa_der(obj.expr);
            bool = reshape(bool, obj.szs{:});
            id = reshape(id, obj.szs{:});
        end
        
        function boolv = isa_reducible(obj)
            if all(isa_reducible(obj.expr))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
            
        function value = evaluate(obj)
            value = reshape(evaluate(obj.expr), obj.szs{:});
        end
        
        function v = forward(obj)
            obj.m_value = reshape(value(obj.expr), obj.szs{:});
            v = obj.m_value;
        end
        
        function ast(obj)
            str = [];
            for k=1:length(obj.szs)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['reshape(expr, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.expr);
            end_child(obj);
            for k=1:(length(obj.szs)-1)
                begin_child(obj);
                ast(obj.szs{k});
                end_child(obj);
            end
            last_child(obj);
            ast(obj.szs{end});
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
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end