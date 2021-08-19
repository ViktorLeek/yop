classdef equality_constraint < yop.node
    % h(t,x,z,u,p) <= 0
    properties
        expr
    end
    methods
        function obj = equality_constraint(expr)
            % h(t,x,z,u,p) == 0
            obj@yop.node();
            obj.expr = expr;
        end
        
        function draw(obj)
            fprintf('equality_constraint(expr)\n');

            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end