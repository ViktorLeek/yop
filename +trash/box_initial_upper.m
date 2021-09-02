classdef box_initial_upper < yop.box_constraint
    % Box constraint - upper bound on the initial value
    
    properties (Constant)
        name = 'box_initial_upper'
    end
    
    methods
        function obj = box_initial_upper(var, bnd)
            obj@yop.box_constraint(var, bnd);
        end
    end
end