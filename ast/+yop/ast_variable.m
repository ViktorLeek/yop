classdef ast_variable < yop.ast_node
    
    methods
        
        function obj = ast_variable(dim)    
            if nargin==0
                obj.dim = [1, 1];
            else
                obj.dim = dim;
            end
        end
        
        function ast(obj)
            fprintf('var\n');
        end
        
    end
end