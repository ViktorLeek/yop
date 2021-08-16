classdef ast_variable < yop.ast_expression
    
    properties
        value
    end
    
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
        
        function value = evaluate(obj)
            value = obj.value;
        end
        
    end
end