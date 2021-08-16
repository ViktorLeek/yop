classdef (InferiorClasses = {?yop.ast_expression}) ast_relation < yop.ast_node & yop.ast_ool
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
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function diff = minus(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function id = uplus(expr)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function neg = uminus(expr)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function prod = times(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function prod = mtimes(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = rdivide(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = ldivide(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = mrdivide(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function frac = mldivide(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function pow = power(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function pow = mpower(lhs, rhs)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = transpose(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = ctranspose(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cat(dim, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = horzcat(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = vertcat(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = reshape(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = repmat(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = heaviside(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = abs(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sqrt(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sin(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cos(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = tan(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = asin(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = acos(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = atan(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sinh(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cosh(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = tanh(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = asinh(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = acosh(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = atanh(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = exp(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = log(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = log10(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = floor(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = ceil(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = erf(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = erfinv(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sign(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = mod(a, m)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = atan2(y, x)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = trace(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sum(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = norm(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = min(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = max(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = sumsqr(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = linspace(x1, x2, n)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cross(A, B, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = det(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = inv(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = pinv(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = dot(A, B, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = expm(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = cumsum(obj, varargin)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = der(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function node = alg(obj)
            yop.error('Method not implemented for class "yop.ast_relation".');
        end
        
        function idx = end(obj, k, n)
            idx = builtin('end', ones(size(obj)), k, n);
        end
        
        function ast(obj)
            fprintf([obj.name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            ast(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            ast(obj.rhs);
            end_child(obj);
        end
        
    end
end