classdef ast_timepoint < yop.ast_expression
    properties
        timepoint
        expr
    end
    methods
        function obj = ast_timepoint(timepoint, expr)
            obj.timepoint = timepoint;
            obj.expr = expr;
        end
        
        function ast(obj)
            fprintf('timepoint(timepoint, expr)\n');
            
            begin_child(obj);
            ast(obj.timepoint);
            end_child(obj);
            
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
        
    end
end