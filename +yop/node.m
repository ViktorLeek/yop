classdef node < handle
    % Implements the basic printing behaviour
    
    properties
        dim = [1, 1] % dimensions of 'size'
    end
    
    properties (Constant)
        % reference to handle is constant, but not the value itself.
        stream = yop.stream_state() % named stream to avoid clash with method.
    end
    
    methods
        
        function obj = node()
        end
        
        function sz = size(obj)
            sz = obj.dim;
        end
    end
    
    methods % AST
        
        function sum = plus(lhs, rhs)
            sum = yop.plus(lhs, rhs);
        end
        
        function diff = minus(lhs, rhs)
            diff = yop.minus(lhs, rhs);
        end
        
        function id = uplus(expr)
            id = yop.uplus(expr);
        end
        
        function neg = uminus(expr)
            neg = yop.uminus(expr);
        end
        
        function prod = times(lhs, rhs)
            prod = yop.times(lhs, rhs);
        end
        
        function prod = mtimes(lhs, rhs)
            prod = yop.mtimes(lhs, rhs);
        end
        
        function frac = rdivide(lhs, rhs)
            frac = yop.rdivide(lhs, rhs);
        end
        
        function frac = ldivide(lhs, rhs)
            frac = yop.ldivide(lhs, rhs);
        end
        
        function frac = mrdivide(lhs, rhs)
            frac = yop.mrdivide(lhs, rhs);
        end
        
        function frac = mldivide(lhs, rhs)
            frac = yop.mldivide(lhs, rhs);
        end
        
        function pow = power(lhs, rhs)
            pow = yop.power(lhs, rhs);
        end
        
        function pow = mpower(lhs, rhs)
            pow = yop.mpower(lhs, rhs);
        end
        
        function rel = lt(lhs, rhs)
            rel = yop.lt(lhs, rhs);
        end
        
        function rel = gt(lhs, rhs)
            rel = yop.gt(lhs, rhs);
        end
        
        function rel = le(lhs, rhs)
            rel = yop.le(lhs, rhs);
        end
        
        function rel = ge(lhs, rhs)
            rel = yop.ge(lhs, rhs);
        end
        
        function rel = ne(lhs, rhs)
            rel = yop.ne(lhs, rhs);
        end
        
        function rel = eq(lhs, rhs)
            rel = yop.eq(lhs, rhs);
        end
        
        function varargout = subsref(obj, s)
            % sr = subsref(obj, s)
            % The function is designed to enable two things.
            %   1) subsref into variables and expressions that builds the 
            %      AST. E.g. the code 'expr(1,3)' extracts an element.
            %   2) Time expressions. E.g. 'expr(t==1)', 'expr(t0)'
            %
            % All other cases are refered to the builtin implementation
            
            switch s(1).type                
                    
                case '()'
                    numeric = true;
                    for k=1:length(s.subs)
                        numeric = isnumeric(s.subs{k}) && numeric;
                    end
                    
                    if length(s) == 1  && length(s.subs) == 1 && ...
                            isa(s.subs{1}, 'yop.node')
                        % timed expression
                        
                    elseif length(s) == 1  && numeric
                        varargout{1} = yop.subsref(obj, s);
                        
                    else
                        % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',obj, s);
                        
                    end
                
                otherwise
                    [varargout{1:nargout}] = builtin('subsref', obj, s);
            end
            
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