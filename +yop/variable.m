classdef variable < yop.node
    
    methods
        
        function obj = variable(dim)
            if nargin==0
                obj.dim = [1, 1];
            else
                obj.dim = dim;
            end
        end
        
        function print(obj)
            fprintf('var\n');
        end
        
    end
end