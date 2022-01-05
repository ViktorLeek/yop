classdef ast_expression < yop.ast_node
    % ast_expression
    % The purpose of this class is to enable to right operator and function
    % overloads for expressions. This class is inferior to the
    % yop.ast_relation class in order to avoid structures with unclear
    % semantics, such as '(expr1 <= expr2) + expr3'
    
    properties
        dim = [1, 1] % dimensions of 'size'
        m_ival = false
    end
    
    methods
        
        function obj = ast_expression(ival)
            obj@yop.ast_node();
            if nargin == 1
                obj.m_ival = ival;
            end
        end
        
        function bool = is_ival(obj)
            % Tests if a node is an interval
            % The reason for having this implemented this way is that
            % intervals are implicitly propagated to child nodes.
            bool = obj.m_ival;
        end
        
        function [t0, tf] = get_ival(obj)
            % This analysis is optimistic, which is bad. This does not
            % consider the case that the interval can be in a subexpression
            % that does not reach the final expression. A fix that tests if
            % the found interval reaches the expression is necessary in
            % order for correct analysis.
            t0 = yop.initial_timepoint;
            tf = yop.final_timepoint;
            [tsort, N] = topological_sort(obj);
            for n=1:N
                if isa(tsort{n}, 'yop.ast_timeinterval')
                    % The function returns here, so having more than one
                    % interval in an expression is undefined behaviour.
                    t0 = tsort{n}.t0;
                    tf = tsort{n}.tf;
                    return
                end
            end
        end
        
        function sz = size(obj, varargin)
            sz = size(ones(obj.dim), varargin{:});
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
            sum = yop.ast_plus(lhs, rhs);
        end
        
        function diff = minus(lhs, rhs)
            diff = yop.ast_minus(lhs, rhs);
        end
        
        function id = uplus(expr)
            id = yop.ast_uplus(expr);
        end
        
        function neg = uminus(expr)
            neg = yop.ast_uminus(expr);
        end
        
        function prod = times(lhs, rhs)
            prod = yop.ast_times(lhs, rhs);
        end
        
        function prod = mtimes(lhs, rhs)
            prod = yop.ast_mtimes(lhs, rhs);
        end
        
        function frac = rdivide(lhs, rhs)
            frac = yop.ast_rdivide(lhs, rhs);
        end
        
        function frac = ldivide(lhs, rhs)
            frac = yop.ast_ldivide(lhs, rhs);
        end
        
        function frac = mrdivide(lhs, rhs)
            frac = yop.ast_mrdivide(lhs, rhs);
        end
        
        function frac = mldivide(lhs, rhs)
            frac = yop.ast_mldivide(lhs, rhs);
        end
        
        function pow = power(lhs, rhs)
            pow = yop.ast_power(lhs, rhs);
        end
        
        function pow = mpower(lhs, rhs)
            pow = yop.ast_mpower(lhs, rhs);
        end
        
        function node = transpose(obj)
            node = yop.ast_transpose(obj);
        end
        
        function node = ctranspose(obj)
            node = yop.ast_ctranspose(obj);
        end
        
        function node = cat(dim, varargin)
            node = yop.ast_cat(dim, varargin{:});
        end
        
        function node = horzcat(obj, varargin)
            node = yop.ast_horzcat(obj, varargin{:});
        end
        
        function node = vertcat(obj, varargin)
            node = yop.ast_vertcat(obj, varargin{:});
        end
        
        function node = reshape(obj, varargin)
            node = yop.ast_reshape(obj, varargin{:});
        end
        
        function node = repmat(obj, varargin)
            node = yop.ast_repmat(obj, varargin{:});
        end
        
        function node = heaviside(obj)
            node = yop.ast_heaviside(obj);
        end
        
        function node = abs(obj)
            node = yop.ast_abs(obj);
        end
        
        function node = sqrt(obj)
            node = yop.ast_sqrt(obj);
        end
        
        function node = sin(obj)
            node = yop.ast_sin(obj);
        end
        
        function node = cos(obj)
            node = yop.ast_cos(obj);
        end
        
        function node = tan(obj)
            node = yop.ast_tan(obj);
        end
        
        function node = asin(obj)
            node = yop.ast_asin(obj);
        end
        
        function node = acos(obj)
            node = yop.ast_acos(obj);
        end
        
        function node = atan(obj)
            node = yop.ast_atan(obj);
        end
        
        function node = sinh(obj)
            node = yop.ast_sinh(obj);
        end
        
        function node = cosh(obj)
            node = yop.ast_cosh(obj);
        end
        
        function node = tanh(obj)
            node = yop.ast_tanh(obj);
        end
        
        function node = asinh(obj)
            node = yop.ast_asinh(obj);
        end
        
        function node = acosh(obj)
            node = yop.ast_acosh(obj);
        end
        
        function node = atanh(obj)
            node = yop.ast_atanh(obj);
        end
        
        function node = exp(obj)
            node = yop.ast_exp(obj);
        end
        
        function node = log(obj)
            node = yop.ast_log(obj);
        end
        
        function node = log10(obj)
            node = yop.ast_log10(obj);
        end
        
        function node = floor(obj)
            node = yop.ast_floor(obj);
        end
        
        function node = ceil(obj)
            node = yop.ast_ceil(obj);
        end
        
        function node = erf(obj)
            node = yop.ast_erf(obj);
        end
        
        function node = erfinv(obj)
            node = yop.ast_erfinv(obj);
        end
        
        function node = sign(obj)
            node = yop.ast_sign(obj);
        end
        
        function node = mod(a, m)
            node = yop.ast_mod(a, m);
        end
        
        function node = atan2(y, x)
            node = yop.ast_atan2(y, x);
        end
        
        function node = trace(obj)
            node = yop.ast_trace(obj);
        end
        
        function node = sum(obj, varargin)
            node = yop.ast_sum(obj, varargin{:});
        end
        
        function node = norm(obj, varargin)
            node = yop.ast_norm(obj, varargin{:});
        end
        
        function node = min(obj, varargin)
            node = yop.ast_min(obj, varargin{:});
        end
        
        function node = max(obj, varargin)
            node = yop.ast_max(obj, varargin{:});
        end
        
        function node = sumsqr(obj)
            node = yop.ast_sumsqr(obj);
        end
        
        function node = linspace(x1, x2, n)
            node = yop.ast_linspace(x1, x2, n);
        end
        
        function node = cross(A, B, varargin)
            node = yop.ast_cross(A, B, varargin{:});
        end
        
        function node = det(obj)
            node = yop.ast_det(obj);
        end
        
        function node = inv(obj)
            node = yop.ast_inv(obj);
        end
        
        function node = pinv(obj)
            node = yop.ast_pinv(obj);
        end
        
        function node = dot(A, B, varargin)
            node = yop.ast_dot(A, B, varargin{:});
        end
        
        function node = expm(obj)
            node = yop.ast_expm(obj);
        end
        
        function node = cumsum(obj, varargin)
            node = yop.ast_cumsum(obj, varargin{:});
        end
        
        function node = der(obj)
            node = yop.ast_der(obj);
        end
        
        function node = alg(obj)
            node = yop.ast_alg(obj);
        end
        
        function node = int(obj)
            node = yop.ast_int(obj);
        end
        
        function idx = end(obj, k, n)
            idx = builtin('end', ones(size(obj)), k, n);
        end
        
        
        function varargout = subsref(obj, s)
            % sr = subsref(obj, s)
            % The function is designed to enable two things.
            %   1) subsref into variables and expressions that builds the
            %      AST. E.g. the code 'expr(1,3)' extracts an element.
            %   2) Time expressions. E.g. 'expr(t==1)', 'expr(t0)'
            %
            % All other cases are refered to the builtin implementation.
            % This means that neasted cases do not work. Should be improved
            % later.
            
            % in the future: matlab.mixin.indexing (package)
            
            % Ugly fix to enable expr.at(t==tp)
            if length(s) > 1 && s(1).type == "()"
                indexing = true;
                for k=1:length(s(1).subs)
                    indexing = (isnumeric(s(1).subs{k}) || ischar(s(1).subs{k}) || islogical(s(1).subs{k}) ) ...
                        && indexing;
                end
                if indexing
                    varargout{1} = subsref(subsref(obj, s(1)), s(2:end));
                    return;
                end
                % if it is not indexing, then we simply continue
            end
            
            if length(s) > 1 || s(1).type ~= "()"
                % This needs to be revised because it only works with
                % objects of size [1, 1]. Must overload 
                % numArgumentsFromSubscript
                [varargout{1:nargout}] = builtin('subsref',obj, s);
                return;
            end
            % From this point length(s) == 1 && s.type == "()"

            % Determine if subs contains only numeric elements or colon.
            % subs: {[s1], [s2], ..., [sN]}
            numeric = true;
            for k=1:length(s.subs)
                numeric = (isnumeric(s.subs{k}) || ischar(s.subs{k}) || islogical(s.subs{k}) ) ...
                    && numeric;
            end
            
            if numeric
                % Case: node(numeric_subs)
                varargout{1} = yop.ast_subsref(obj, s);
                
            elseif length(s.subs)==1 && isa(s.subs{1}, 'yop.ast_node')
                % timed expression: obj(t==4)
                % Error handling in ast_timepoint class, as it easier to
                % intuitively understand where that logic is placed if it
                % is inside the constructor of that class.
                varargout{1} = yop.ast_expression.timed_expression( ...
                    s.subs{1}, obj);
                
            else
                % Use built-in for any other expression
                [varargout{1:nargout}] = builtin('subsref',obj, s);
                
            end

        end
        
        function obj = subsasgn(obj, s, varargin)
            
            if length(s) > 1 || s(1).type ~= "()"
                obj = builtin('subsasgn', obj, s, varargin{:});
                return;
            end
            
            numeric = true;
            for k=1:length(s.subs)
                numeric = (isnumeric(s.subs{k}) || ischar(s.subs{k}) ) ...
                    && numeric;
            end
            
            if numeric
                % here it is assumed that varargin{:} is simply a third
                % argument 'b' (see Matlab docs).
                obj = yop.ast_subsasgn(obj, s, varargin{:});
                
            else
                % Use built-in for any other expression
                obj = builtin('subsasgn', obj, s, varargin{:});
                
            end
            
        end  
        
    end
    
    methods (Static)
        function ast = timed_expression(texpr, expr)
            
            ist0 = isa_independent0(texpr); % expr(t0)
            istf = isa_independentf(texpr); % expr(tf)
            ist  = isa_independent(texpr);  % expr(t)
            
            isrel = isa(texpr, 'yop.ast_relation');
            ssr = yop.to_ssr(texpr);
            n = length(ssr);
            
            if isrel && n==1
                lhs1 = ssr{1}.lhs;
                rhs1 = ssr{1}.rhs;
                a1_is_eq = isa(ssr{1}, 'yop.ast_eq');
                a1_is_le = isa(ssr{1}, 'yop.ast_lt') || isa(ssr{1}, 'yop.ast_le');
                a1_is_ge = isa(ssr{1}, 'yop.ast_gt') || isa(ssr{1}, 'yop.ast_ge');
                a1_lhs_is_t = isa_independent(lhs1);
                a1_rhs_is_t = isa_independent(rhs1);
                a1_lhs_is_t0 = isa_independent0(lhs1);
                a1_rhs_is_t0 = isa_independent0(rhs1);
                a1_lhs_is_tf = isa_independentf(lhs1);
                a1_rhs_is_tf = isa_independentf(rhs1);
                a1_lhs_is_num = isa_numeric(lhs1);
                a1_rhs_is_num = isa_numeric(rhs1);
                
                
            elseif isrel && n==2
                lhs1 = ssr{1}.lhs;
                rhs1 = ssr{1}.rhs;
                a1_is_le = isa(ssr{1}, 'yop.ast_lt') || isa(ssr{1}, 'yop.ast_le');
                a1_is_ge = isa(ssr{1}, 'yop.ast_gt') || isa(ssr{1}, 'yop.ast_ge');
                a1_lhs_is_t = isa_independent(lhs1);
                a1_rhs_is_t0 = isa_independent0(rhs1);
                a1_rhs_is_tf = isa_independentf(rhs1);
                a1_rhs_is_num = isa_numeric(rhs1);
                
                lhs2 = ssr{2}.lhs;
                rhs2 = ssr{2}.rhs;
                a2_is_le = isa(ssr{2}, 'yop.ast_lt') || isa(ssr{2}, 'yop.ast_le');
                a2_is_ge = isa(ssr{2}, 'yop.ast_gt') || isa(ssr{2}, 'yop.ast_ge');
                a2_rhs_is_t = isa_independent(rhs2);
                a2_lhs_is_t0 = isa_independent0(lhs2);
                a2_lhs_is_tf = isa_independentf(lhs2);
                a2_lhs_is_num = isa_numeric(lhs2);
                
            end
            
            % Notice that the comparison to n (which is a value that always
            % exist here) is done first. Since Matlab applies short 
            % circuiting none of the other logical expressions are 
            % evaluated if the first one failes, so it does not matter if
            % several of the others may not have a definition.
            if ist0 
                % expr(t0)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif istf
                % expr(tf)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif ist
                % expr(t)
                ast = expr;
                
            elseif n==1 && a1_lhs_is_t && a1_is_eq && a1_rhs_is_t0
                % expr(t == t0)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_t && a1_is_eq && a1_rhs_is_tf
                % expr(t == tf)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_t && a1_is_eq && a1_rhs_is_num
                % expr(t == num)
                ast = yop.ast_timepoint(yop.prop_num(rhs1), expr);
                
            elseif n==1 && a1_lhs_is_t0 && a1_is_eq && a1_rhs_is_t
                % expr(t0 == t)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_tf && a1_is_eq && a1_rhs_is_t
                % expr(tf == t)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_num && a1_is_eq && a1_rhs_is_t
                % expr(num == t)
                ast = yop.ast_timepoint(yop.prop_num(lhs1), expr);
                
            elseif n==1 && a1_lhs_is_t && a1_is_le && a1_rhs_is_t0
                % expr(t <= t0)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_t && a1_is_le && a1_rhs_is_tf
                % expr(t <= tf)
                ast = expr;
                
            elseif n==1 && a1_lhs_is_t && a1_is_le && a1_rhs_is_num
                % expr(t <= num)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    yop.prop_num(rhs1), ...
                    expr);
                
            elseif n==1 && a1_lhs_is_t0 && a1_is_le && a1_rhs_is_t
                % expr(t0 <= t)
                ast = expr;
                
            elseif n==1 && a1_lhs_is_tf && a1_is_le && a1_rhs_is_t
                % expr(tf <= t)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_num && a1_is_le && a1_rhs_is_t
                % expr(num <= t)
                ast = yop.ast_timeinterval( ...
                    yop.prop_num(lhs1), ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==1 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_t0
                % expr(t >= t0)
                ast = expr;
                
            elseif n==1 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_tf
                % expr(t >= tf)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_num
                % expr(t >= num)
                ast = yop.ast_timeinterval( ...
                    yop.prop_num(rhs1), ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==1 && a1_lhs_is_t0 && a1_is_ge && a1_rhs_is_t
                % expr(t0 >= t)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_lhs_is_tf && a1_is_ge && a1_rhs_is_t
                % expr(tf >= t)
                ast = expr;
                
            elseif n==1 && a1_lhs_is_num && a1_is_ge && a1_rhs_is_t
                % expr(num >= t)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    yop.prop_num(lhs1), ...
                    expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_le && a1_rhs_is_t0 && a2_lhs_is_t0 && a2_is_le && a2_rhs_is_t
                % expr(t0 <= t <= t0) -> expr(t <= t0), expr(t0 <= t)
                % Operator binding and dfs gives the two relations
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_le && a1_rhs_is_tf && a2_lhs_is_t0 && a2_is_le && a2_rhs_is_t
                % expr(t0 <= t <= tf) -> expr(t <= tf), expr(t0 <= t)
                ast = expr;
                
            elseif n==2 && a1_lhs_is_t && a1_is_le && a1_rhs_is_tf && a2_lhs_is_tf && a2_is_le && a2_rhs_is_t
                % expr(tf <= t <= tf) -> expr(t <= tf), expr(tf <= t)
                % Operator binding and dfs gives the two relations
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_le && a1_rhs_is_tf && a2_lhs_is_num && a2_is_le && a2_rhs_is_t
                % expr(num <= t <= tf) -> expr(t <= tf), expr(num <= t)
                ast = yop.ast_timeinterval( ...
                    yop.prop_num(lhs2), ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_le && a1_rhs_is_num && a2_lhs_is_t0 && a2_is_le && a2_rhs_is_t
                % expr(t0 <= t <= num) -> expr(t <= num), expr(t0 <= t)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    yop.prop_num(rhs1), ...
                    expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_le && a1_rhs_is_num && a2_lhs_is_num && a2_is_le && a2_rhs_is_t
                % expr(num <= t <= num) -> expr(t <= num), expr(num <= t)
                ast = yop.ast_timeinterval( ...
                    yop.prop_num(lhs2), ...
                    yop.prop_num(rhs1), ...
                    expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_t0 && a2_lhs_is_t0 && a2_is_ge && a2_rhs_is_t
                % expr(t0 >= t >= t0) -> expr(t >= t0), expr(t0 >= t)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_t0 && a2_lhs_is_tf && a2_is_ge && a2_rhs_is_t
                % expr(tf >= t >= t0) -> expr(t >= t0), expr(tf >= t)
                ast = expr;
                
            elseif n==2 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_t0 && a2_lhs_is_num && a2_is_ge && a2_rhs_is_t
                % expr(num >= t >= t0) -> expr(t >= t0), expr(num >= t)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    yop.prop_num(lhs2), ...
                    expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_tf && a2_lhs_is_tf && a2_is_ge && a2_rhs_is_t
                % expr(tf >= t >= tf) -> expr(t >= tf), expr(tf >= t)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_num && a2_lhs_is_tf && a2_is_ge && a2_rhs_is_t
                % expr(tf >= t >= num) -> expr(t >= num), expr(tf >= t)
                ast = yop.ast_timeinterval( ...
                    yop.prop_num(rhs1), ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==2 && a1_lhs_is_t && a1_is_ge && a1_rhs_is_num && a2_lhs_is_num && a2_is_ge && a2_rhs_is_t
                % expr(num >= t >= num) -> expr(t >= num), expr(num >= t)
                ast = yop.ast_timeinterval( ...
                    yop.prop_num(rhs1), ...
                    yop.prop_num(lhs2), ...
                    expr);
                
            else
                error(yop.error.failed_to_parse_timed_expression());
                
            end
            
        end
    end
end