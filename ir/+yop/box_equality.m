classdef box_equality < yop.box_constraint
    
    properties (Constant)
        name = 'box_eq'
    end
    
    methods
        function obj = box_equality(var, bnd)
            obj@yop.box_constraint(var, bnd);
        end
    end
end