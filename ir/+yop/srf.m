classdef srf < yop.node
    % srf - single relation form
    properties
        lhs
        rhs
    end
    methods
        function obj = srf(lhs, rhs)
            % A relation on the form 'lhs' 'rel' 'rhs'. E.g.
            %  exp1 <= exp2
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
        
        function draw(obj)
            fprintf([obj.name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            draw(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            draw(obj.rhs);
            end_child(obj);
        end
    end
end