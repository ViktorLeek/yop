classdef OOL < handle
    % Operator overloading. Inherited by the classes that make up the AST.
    methods
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
        
    end
end