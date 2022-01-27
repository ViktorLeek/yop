classdef ast_expression < yop.ast_node
    % ast_expression
    % The purpose of this class is to enable to right operator and function
    % overloads for expressions. This class is inferior to the
    % yop.ast_relation class in order to avoid structures with confusing
    % semantics, such as '(expr1 <= expr2) + expr3'
    
    properties
        m_numval
        m_t0
        m_tf
        m_der
        m_reducible
        m_type
        m_typeid
    end
    
    methods
        
        function obj = ast_expression(value, numval, t0, tf, isder, isreducible, type, typeid)
            obj@yop.ast_node(value);
            obj.m_numval    = numval;
            obj.m_t0        = t0;
            obj.m_tf        = tf;
            obj.m_der       = isder;
            obj.m_reducible = isreducible;
            obj.m_type      = type;
            obj.m_typeid    = typeid;
        end
        
        function val = numval(obj)
            val = obj.m_numval;
        end
        
        function bool = isa_numeric(obj)
            bool = ~isnan(obj.m_numval);
        end
        
        function t0 = get_t0(obj)
            t0 = obj.m_t0;
        end
        
        function tf = get_tf(obj)
            tf = obj.m_tf;
        end
        
        function bool = isa_ival(obj)
            bool = ~isinf(obj.m_t0) | ~isinf(obj.m_tf);
        end
        
        function [t0, tf] = get_ival(obj)
            t0 = obj.m_t0;
            tf = obj.m_tf;
        end
        
        function bool = isa_timepoint(obj)
            bool = obj.m_t0 == obj.m_tf;
        end
        
        function id = get_der(obj)
            id = obj.m_der;
        end
        
        function bool = isa_der(obj)
            bool = obj.m_der > 0;
        end
        
        function boolv = isa_reducible(obj)
            boolv = obj.m_reducible;
        end
        
        function [type, id] = Type(obj)
            type = obj.m_type;
            id   = obj.m_typeid;
        end
        
        function bool = isa_variable(obj)
            bool = ...
                obj.m_type >= yop.var_type.variables_start & ...
                obj.m_type <= yop.var_type.variables_stop;
        end
        
        function sz = size(obj, varargin)
            sz = size(obj.m_numval, varargin{:});
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
        
        function node = int(obj)
            node = yop.ast_int(obj);
        end
        
        function node = if_else(varargin)
            node = yop.ast_if_else(varargin{:});
        end
        
        function idx = end(obj, k, n)
            idx = builtin('end', ones(size(obj)), k, n);
        end
        
        function node = at(expression, timepoint)
            % Alternative syntax for evaluating expression at a timepoint
            %node = yop.ast_timepoint(timepoint, expression);
            node=yop.ast_expression.timed_expression(timepoint, expression);
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
            
            type = Type(texpr);
            ist  = type == yop.var_type.time;  % expr(t)
            ist0 = type == yop.var_type.time0; % expr(t0)
            istf = type == yop.var_type.timef; % expr(tf)
            
            isrel = isa(texpr, 'yop.ast_relation');
            ssr = yop.to_ssr(texpr);
            n = length(ssr);
            
            if isrel && n==1
                lhs1 = ssr{1}.m_lhs;
                rhs1 = ssr{1}.m_rhs;
                
                a1_eq = isa(ssr{1}, 'yop.ast_eq');
                a1_le = isa(ssr{1}, 'yop.ast_lt') || isa(ssr{1}, 'yop.ast_le');
                a1_ge = isa(ssr{1}, 'yop.ast_gt') || isa(ssr{1}, 'yop.ast_ge');
                
                a1_isnum_lhs = isa_numeric(lhs1);
                a1_isnum_rhs = isa_numeric(rhs1);
                a1_num_lhs = numval(lhs1);
                a1_num_rhs = numval(rhs1);
                
                a1_type_lhs = Type(lhs1);
                a1_type_rhs = Type(rhs1);
                
                a1_time_lhs  = a1_type_lhs == yop.var_type.time;
                a1_time0_lhs = a1_type_lhs == yop.var_type.time0;
                a1_timef_lhs = a1_type_lhs == yop.var_type.timef;
                
                a1_time_rhs  = a1_type_rhs == yop.var_type.time;
                a1_time0_rhs = a1_type_rhs == yop.var_type.time0;
                a1_timef_rhs = a1_type_rhs == yop.var_type.timef;
                
            elseif isrel && n==2
                lhs1 = ssr{1}.m_lhs;
                rhs1 = ssr{1}.m_rhs;
                
                a1_le = isa(ssr{1}, 'yop.ast_lt') || isa(ssr{1}, 'yop.ast_le');
                a1_ge = isa(ssr{1}, 'yop.ast_gt') || isa(ssr{1}, 'yop.ast_ge');
                
                a1_isnum_rhs = isa_numeric(rhs1);
                a1_num_rhs = numval(rhs1);
                
                a1_type_lhs = Type(lhs1);
                a1_type_rhs = Type(rhs1);
                
                a1_time_lhs  = a1_type_lhs == yop.var_type.time;
                a1_time0_rhs = a1_type_rhs == yop.var_type.time0;
                a1_timef_rhs = a1_type_rhs == yop.var_type.timef;
                
                lhs2 = ssr{2}.m_lhs;
                rhs2 = ssr{2}.m_rhs;
                
                a2_le = isa(ssr{2}, 'yop.ast_lt') || isa(ssr{2}, 'yop.ast_le');
                a2_ge = isa(ssr{2}, 'yop.ast_gt') || isa(ssr{2}, 'yop.ast_ge');
                
                a2_type_lhs = Type(lhs2);
                a2_type_rhs = Type(rhs2);
                
                a2_time0_lhs = a2_type_lhs == yop.var_type.time0;
                a2_timef_lhs = a2_type_lhs == yop.var_type.timef;
                a2_isnum_lhs = isa_numeric(lhs2);
                a2_num_lhs = numval(lhs2);
                
                a2_time_rhs = a2_type_rhs == yop.var_type.time;
                
            end
            
            % Notice that the comparison to n (which is a value that is
            % defined here) is done first. Since Matlab applies short 
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
                
            elseif n==1 && a1_time_lhs && a1_eq && a1_time0_rhs
                % expr(t == t0)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_time_lhs && a1_eq && a1_timef_rhs
                % expr(t == tf)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_time_lhs && a1_eq && a1_isnum_rhs
                % expr(t == num)
                ast = yop.ast_timepoint(a1_num_rhs, expr);
                
            elseif n==1 && a1_time0_lhs && a1_eq && a1_time_rhs
                % expr(t0 == t)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_timef_lhs && a1_eq && a1_time_rhs
                % expr(tf == t)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_isnum_lhs && a1_eq && a1_time_rhs
                % expr(num == t)
                ast = yop.ast_timepoint(a1_num_lhs, expr);
                
            elseif n==1 && a1_time_lhs && a1_le && a1_time0_rhs
                % expr(t <= t0)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_time_lhs && a1_le && a1_timef_rhs
                % expr(t <= tf)
                ast = expr;
                
            elseif n==1 && a1_time_lhs && a1_le && a1_isnum_rhs
                % expr(t <= num)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    a1_num_rhs, ...
                    expr);
                
            elseif n==1 && a1_time0_lhs && a1_le && a1_time_rhs
                % expr(t0 <= t)
                ast = expr;
                
            elseif n==1 && a1_timef_lhs && a1_le && a1_time_rhs
                % expr(tf <= t)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_isnum_lhs && a1_le && a1_time_rhs
                % expr(num <= t)
                ast = yop.ast_timeinterval( ...
                    a1_num_lhs, ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==1 && a1_time_lhs && a1_ge && a1_time0_rhs
                % expr(t >= t0)
                ast = expr;
                
            elseif n==1 && a1_time_lhs && a1_ge && a1_timef_rhs
                % expr(t >= tf)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==1 && a1_time_lhs && a1_ge && a1_isnum_rhs
                % expr(t >= num)
                ast = yop.ast_timeinterval( ...
                    a1_num_rhs, ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==1 && a1_time0_lhs && a1_ge && a1_time_rhs
                % expr(t0 >= t)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==1 && a1_timef_lhs && a1_ge && a1_time_rhs
                % expr(tf >= t)
                ast = expr;
                
            elseif n==1 && a1_isnum_lhs && a1_ge && a1_time_rhs
                % expr(num >= t)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    a1_num_lhs, ...
                    expr);
                
            elseif n==2 && a1_time_lhs && a1_le && a1_time0_rhs && a2_time0_lhs && a2_le && a2_time_rhs
                % expr(t0 <= t <= t0) -> expr(t <= t0), expr(t0 <= t)
                % Operator binding and dfs gives the two relations
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==2 && a1_time_lhs && a1_le && a1_timef_rhs && a2_time0_lhs && a2_le && a2_time_rhs
                % expr(t0 <= t <= tf) -> expr(t <= tf), expr(t0 <= t)
                ast = expr;
                
            elseif n==2 && a1_time_lhs && a1_le && a1_timef_rhs && a2_timef_lhs && a2_le && a2_time_rhs
                % expr(tf <= t <= tf) -> expr(t <= tf), expr(tf <= t)
                % Operator binding and dfs gives the two relations
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==2 && a1_time_lhs && a1_le && a1_timef_rhs && a2_isnum_lhs && a2_le && a2_time_rhs
                % expr(num <= t <= tf) -> expr(t <= tf), expr(num <= t)
                ast = yop.ast_timeinterval( ...
                    a2_num_lhs, ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==2 && a1_time_lhs && a1_le && a1_isnum_rhs && a2_time0_lhs && a2_le && a2_time_rhs
                % expr(t0 <= t <= num) -> expr(t <= num), expr(t0 <= t)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    a1_num_rhs, ...
                    expr);
                
            elseif n==2 && a1_time_lhs && a1_le && a1_isnum_rhs && a2_isnum_lhs && a2_le && a2_time_rhs
                % expr(num <= t <= num) -> expr(t <= num), expr(num <= t)
                ast = yop.ast_timeinterval( ...
                    a2_num_lhs, ...
                    a1_num_rhs, ...
                    expr);
                
            elseif n==2 && a1_time_lhs && a1_ge && a1_time0_rhs && a2_time0_lhs && a2_ge && a2_time_rhs
                % expr(t0 >= t >= t0) -> expr(t >= t0), expr(t0 >= t)
                ast = yop.ast_timepoint(yop.initial_timepoint, expr);
                
            elseif n==2 && a1_time_lhs && a1_ge && a1_time0_rhs && a2_timef_lhs && a2_ge && a2_time_rhs
                % expr(tf >= t >= t0) -> expr(t >= t0), expr(tf >= t)
                ast = expr;
                
            elseif n==2 && a1_time_lhs && a1_ge && a1_time0_rhs && a2_isnum_lhs && a2_ge && a2_time_rhs
                % expr(num >= t >= t0) -> expr(t >= t0), expr(num >= t)
                ast = yop.ast_timeinterval( ...
                    yop.initial_timepoint, ...
                    a2_num_lhs, ...
                    expr);
                
            elseif n==2 && a1_time_lhs && a1_ge && a1_timef_rhs && a2_timef_lhs && a2_ge && a2_time_rhs
                % expr(tf >= t >= tf) -> expr(t >= tf), expr(tf >= t)
                ast = yop.ast_timepoint(yop.final_timepoint, expr);
                
            elseif n==2 && a1_time_lhs && a1_ge && a1_isnum_rhs && a2_timef_lhs && a2_ge && a2_time_rhs
                % expr(tf >= t >= num) -> expr(t >= num), expr(tf >= t)
                ast = yop.ast_timeinterval( ...
                    a1_num_rhs, ...
                    yop.final_timepoint, ...
                    expr);
                
            elseif n==2 && a1_time_lhs && a1_ge && a1_isnum_rhs && a2_isnum_lhs && a2_ge && a2_time_rhs
                % expr(num >= t >= num) -> expr(t >= num), expr(num >= t)
                ast = yop.ast_timeinterval( ...
                    a1_num_rhs, ...
                    a2_num_lhs, ...
                    expr);
                
            else
                error(yop.error.failed_to_parse_timed_expression());
                
            end
            
        end
    end
end