classdef srf < yop.ast_node
    % srf - single relation form
    % It derives from ast_node in order to be draw()-able, otherwise
    % intended to be separate from the ast-representation.
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