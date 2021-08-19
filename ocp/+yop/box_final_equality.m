classdef box_final_equality < yop.box_constraint
    % Box constraint - upper bound on the initial value
    
    properties (Constant)
        name = 'box_final_equality'
    end
    
    methods
        function obj = box_final_equality(var, bnd)
            obj@yop.box_constraint(var, bnd);
        end
    end
end