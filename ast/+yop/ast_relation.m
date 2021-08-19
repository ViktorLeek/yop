classdef (InferiorClasses = {?yop.ast_expression, ?yop.ast_variable}) ast_relation < yop.node & yop.ast_ool
    % yop.ast_relation
    % yop.ast_expression is inferior in order to make it impossible to mix
    % expressions and relations in illegal ways, for instance:
    % "(expr <= 4) + 2"
    properties
        lhs
        rhs
    end
    
    methods
        function obj = ast_relation(lhs, rhs)
            obj@yop.node();
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
        
        function rel = lt(lhs, rhs)
            rel = yop.ast_lt(lhs, rhs);
        end
        
        function rel = gt(lhs, rhs)
            rel = yop.ast_gt(lhs, rhs);
        end
        
        function rel = le(lhs, rhs)
            rel = yop.ast_le(lhs, rhs);
        end
        
        function rel = ge(lhs, rhs)
            rel = yop.ast_ge(lhs, rhs);
        end
        
        function rel = ne(lhs, rhs)
            rel = yop.ast_ne(lhs, rhs);
        end
        
        function rel = eq(lhs, rhs)
            rel = yop.ast_eq(lhs, rhs);
        end
        
        function sum = plus(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function diff = minus(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function id = uplus(expr)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function neg = uminus(expr)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function prod = times(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function prod = mtimes(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = rdivide(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = ldivide(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = mrdivide(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = mldivide(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function pow = power(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function pow = mpower(lhs, rhs)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = transpose(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = ctranspose(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cat(dim, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = horzcat(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = vertcat(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = reshape(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = repmat(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = heaviside(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = abs(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sqrt(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sin(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cos(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = tan(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = asin(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = acos(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = atan(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sinh(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cosh(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = tanh(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = asinh(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = acosh(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = atanh(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = exp(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = log(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = log10(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = floor(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = ceil(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = erf(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = erfinv(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sign(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = mod(a, m)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = atan2(y, x)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = trace(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sum(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = norm(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = min(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = max(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sumsqr(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = linspace(x1, x2, n)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cross(A, B, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = det(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = inv(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = pinv(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = dot(A, B, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = expm(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cumsum(obj, varargin)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = der(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function node = alg(obj)
            error('[yop] Error: Method not implemented for class "yop.ast_relation".');
        end
        
        function idx = end(obj, k, n)
            idx = builtin('end', ones(size(obj)), k, n);
        end
        
        function rels = get_relations(obj)
            l = get_relations(obj.lhs);
            r = get_relations(obj.rhs);
            rels = {obj, l{:}, r{:}}; % [{obj}, l(:)', r(:)'];
            rels = rels(~cellfun('isempty', rels));
        end
        
        function r = rmost(obj)
            r = rmost(obj.rhs);
        end
        
        function l = lmost(obj)
            l = lmost(obj.lhs);
        end
        
        function draw(obj)
            fprintf([obj.name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            draw(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            draw(obj.rhs);
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
                topological_sort(obj.lhs, topsort, visited, n_elem);
            
            [topsort, visited, n_elem] = ...
                topological_sort(obj.rhs, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end