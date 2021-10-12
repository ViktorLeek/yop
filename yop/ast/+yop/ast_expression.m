classdef ast_expression < yop.node & yop.ast_ool
    % ast_expression
    % The purpose of this class is to enable to right operator and function
    % overloads for expressions. This class is inferior to the
    % yop.ast_relation class in order to avoid structures with unclear
    % semantics, such as '(expr1 <= expr2) + expr3'
    
    properties
        dim = [1, 1] % dimensions of 'size'
    end
    
    methods
        function [bool, t0, tf] = isa_timeinterval(obj)
            t0 = yop.initial_timepoint;
            tf = yop.final_timepoint;
            bool = false;
            [tsort, N] = topological_sort(obj);
            for n=1:N
                if isa(tsort{n}, 'yop.ast_timeinterval')
                    % The function returns here, so having more than one
                    % interval in an expression is undefined behaviour.
                    bool = true;
                    t0 = tsort{n}.t0;
                    tf = tsort{n}.tf;
                    return
                end
            end
        end
    end
    
    methods
        
        function obj = ast_expression()
            obj@yop.node();
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
                
            elseif length(s.subs)==1 && isa(s.subs{1}, 'yop.node')
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
            
            if isa(texpr, 'yop.ast_independent_initial') || ...
                    isa(texpr, 'yop.ast_independent_final')
                ast = yop.ast_timepoint(texpr, expr);
                return;
            end
            
            if isa(texpr, 'yop.ast_independent')
                ast = expr; % Nothing to be done, expr(t) == expr
                return;
            end
            
            % Slimy, but filling.
            if isa(texpr, 'yop.ast_relation')
                srf = yop.ocp.to_srf({texpr});
                switch length(srf)
                    case 1
                        switch class(srf{1})
                            case 'yop.ast_eq'
                                ast = yop.ast_timepoint(texpr, expr);
                                
                            case {'yop.ast_lt', 'yop.ast_le'}
                                if isa(srf{1}.lhs, 'yop.ast_independent')
                                    % t < value
                                    ast = yop.ast_timeinterval( ...
                                        yop.initial_timepoint, ...
                                        srf{1}.rhs, expr);
                                    
                                else % value < t
                                    ast = yop.ast_timeinterval( ...
                                        srf{1}.lhs, ...
                                        yop.final_timepoint,...
                                        expr);
                                end
                            case {'yop.ast_gt', 'yop.ast_ge'}
                                if isa(srf{1}.lhs, 'yop.ast_independent')
                                    % t > value
                                    ast = yop.ast_timeinterval( ...
                                        srf{1}.rhs, ...
                                        yop.final_timepoint,...
                                        expr);
                                    
                                else % value > t
                                    ast = yop.ast_timeinterval( ...
                                        yop.initial_timepoint, ...
                                        srf{1}.lhs, expr);
                                end
                        end
                        
                    case 2
                        t0=[]; tf=[];
                        for k=1:2
                            switch class(srf{k})
                                case {'yop.ast_lt', 'yop.ast_le'}
                                    if isa(srf{k}.lhs, ...
                                            'yop.ast_independent')
                                        % t < value
                                        tf = srf{k}.rhs;
                                    else % value < t. Should be an elseif to test rhs for independent
                                        t0 = srf{k}.lhs;
                                    end
                                case {'yop.ast_gt', 'yop.ast_ge'}
                                    if isa(srf{k}.lhs, ...
                                            'yop.ast_independent')
                                        % t > value
                                        t0 = srf{k}.rhs;
                                        
                                    else % value > t
                                        tf = srf{k}.lhs;
                                    end
                                    
                                otherwise
                                    error(yop.msg.ival_relation_error(...
                                        class(srf{k})));
                            end
                        end
                        if isempty(t0) || isempty(tf)
                            error(yop.msg.ambig_ival);
                        end
                        ast = yop.ast_timeinterval(t0, tf, expr);
                        
                    otherwise
                        error(yop.msg.timed_expr_relation_overflow);
                end
                
            else
                error(yop.msg.illegal_timepoint);
            end
            
        end
    end
end