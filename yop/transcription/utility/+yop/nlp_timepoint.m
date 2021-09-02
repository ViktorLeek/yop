classdef nlp_timepoint < handle
    properties
        node        
        value = [];
    end
    methods
        function obj = nlp_timepoint(ast_timepoint)
            obj.node = ast_timepoint;                      
        end
        
        function t = tp(obj)
            t = obj.node.timepoint;
        end
        
        function e = expr(obj)
            e = obj.node.expr;
        end
        
        function obj = parameterize(obj)
            for k=1:length(obj)
                obj(k).node.m_value = obj(k).value;
            end
        end
    end
end