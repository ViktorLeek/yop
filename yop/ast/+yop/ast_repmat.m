classdef ast_repmat < yop.ast_expression
    properties
        expr
        args
    end
    methods
        function obj = ast_repmat(expr, varargin)
            obj@yop.ast_expression(is_ival(expr));
            obj.expr = expr;
            obj.args = varargin;
            obj.dim = repmat(ones(size(expr)), varargin{:});
        end
        
        function val = numval(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = numval(obj.args{k});
            end
            val = repmat(numval(obj.expr), tmp{:});
        end
        
        function [type, id] = Type(obj)
            if ~isempty(obj.m_type)
                type = obj.m_type.type;
                id = obj.m_type.id;
                return;
            end
            [type, id] = Type(obj.expr);
            type = repmat(type, obj.args{:});
            id = repmat(id, obj.args{:});
            obj.m_type.type = type;
            obj.m_type.id = id;
        end
        
        function boolv = isa_reducible(obj)
            if all(isa_reducible(obj.expr))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function value = evaluate(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = evaluate(obj.args{k});
            end
            value = repmat(evaluate(obj.expr), tmp{:});
        end
        
        function v = forward(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = value(obj.args{k});
            end
            obj.m_value = repmat(value(obj.expr), tmp{:});
            v = obj.m_value;
        end
        
        function ast(obj)
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['repmat(A, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.expr);
            end_child(obj);
            for k=1:(length(obj.args)-1)
                begin_child(obj);
                ast(obj.args{k});
                end_child(obj);
            end
            last_child(obj);
            ast(obj.args{end});
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
            for k=1:length(obj.args)
                % probably unnecessary, as args are expected to be numerics
                % but could change in the future.
                [topsort, n_elem, visited] = topological_sort( ...
                    obj.args{k}, ...
                    visited, ...
                    topsort, ...
                    n_elem ...
                    );
            end
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end