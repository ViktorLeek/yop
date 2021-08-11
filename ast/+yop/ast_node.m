classdef ast_node < handle
    % Implements the basic printing behaviour
    
    properties
        dim = [1, 1] % dimensions of 'size'
    end
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        stream = yop.stream_state() % named stream to avoid clash with method.
    end
    
    methods
        
        function obj = ast_node()
        end
        
        function sz = size(obj, varargin)
            sz = size(ones(obj.dim), varargin{:});
        end
        
        function ne = numel(obj)
            ne = prod(size(obj));
        end
    end
    
    methods % AST
        
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
            
            if length(s) > 1 || s(1).type ~= "()"
                [varargout{1:nargout}] = builtin('subsref',obj, s);
                return;
            end
            % From this point length(s) == 1 && s.type == "()"

            % Determine if subs contains only numeric elements or colon.
            % subs: {[s1], [s2], ..., [sN]}
            numeric = true;
            for k=1:length(s.subs)
                numeric = (isnumeric(s.subs{k}) || s.subs{k}==":" ) ...
                    && numeric;
            end
            
            if numeric
                % Case: node(numeric_subs)
                varargout{1} = yop.ast_subsref(obj, s);
                
            elseif length(s.subs)==1 && isa(s.subs{1}, 'yop.ast_eq')
                % timed expression: obj(t==4)
                varargout{1} = yop.ast_timepoint(s.subs{1}, obj);
                
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
                numeric = isnumeric(s.subs{k}) && numeric;
            end
            
            if numeric
                % Implement obj(indices) = varargin{:};
                obj = yop.ast_subsasgn(obj, s, varargin{:});
                
            else
                % Use built-in for any other expression
                obj = builtin('subsasgn', obj, s, varargin{:});
                
            end
            
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
        
        function idx = end(obj, k, n)
            idx = builtin('end', ones(size(obj)), k, n);
        end
        
        function node = heaviside(obj)
            node = yop.ast_heaviside(obj);
        end
        
    end
    
    methods % printing
        
        function reset_stream(obj)
            reset(obj.stream);
        end
        
        function indent(obj)
            % Indents the print by implementing the behvaiour:
            % for k=1:obj.stream.indent_level
            %     if obj.stream.branches(k)
            %         fprintf('|');
            %     else
            %         fprintf(' ');
            %     end
            % end
            
            str = repmat(' ', 1, obj.stream.indent_level-1);
            str(obj.stream.branches) = '|';
            fprintf(str);
        end
        
        function indent_more(obj)
            obj.stream.indent_level = obj.stream.indent_level + 2;
        end
        
        function indent_less(obj)
            obj.stream.indent_level = obj.stream.indent_level - 2;
        end
        
        function begin_child(obj)
            indent(obj);
            fprintf('+-');
            obj.stream.branches(obj.stream.indent_level) = true;
            indent_more(obj);
        end
        
        function end_child(obj)
            indent_less(obj)
            obj.stream.branches(obj.stream.indent_level) = false;
        end
        
        function last_child(obj)
            indent(obj);
            fprintf('+-');
            obj.stream.branches(obj.stream.indent_level) = false;
            indent_more(obj);
        end
    end
end