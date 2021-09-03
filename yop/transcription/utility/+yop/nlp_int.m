classdef nlp_int < handle
    properties
        node        
        value = []
        integrand_expr
        integrand_fn
    end
    methods
        function obj = nlp_int(ast_int)
            obj.node = ast_int;                      
        end
        
        function e = integrand(obj)
            e = obj.node.expr;
        end
        
        function obj = parameterize(obj)
            for k=1:length(obj)
                obj(k).node.m_value = obj(k).value;
            end
        end
    end
end