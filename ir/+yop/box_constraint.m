classdef box_constraint < yop.node
    properties
        var
        bnd
    end
    methods
        function obj = box_constraint(var, bnd)
            obj.var = var;
            obj.bnd = bnd;
        end
        
        function draw(obj)
            fprintf([obj.name, '(var, bnd)\n']);
            
            begin_child(obj);
            draw(obj.var);
            end_child(obj);
            
            last_child(obj);
            draw(obj.bnd);
            end_child(obj);
        end
    end
end