classdef box_equality < yop.box_constraint
    
    properties (Constant)
        name = 'box_eq'
    end
    
    methods
        function obj = box_equality(variable, bound)
            obj@yop.box_constraint(variable, bound);
        end
    end
end