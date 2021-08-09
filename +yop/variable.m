classdef variable < handle & OOL
    
    methods
        
        function obj = variable()
        end
        
        function add = plus(lhs, rhs)
            add = yop.add(lhs, rhs);
        end
        
        function sub = minus(lhs, rhs)
            sub = yop.sub(lhs, rhs);
        end
        
        function print(obj)
            fprintf('var\n');
        end
        
    end
end