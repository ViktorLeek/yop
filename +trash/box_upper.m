classdef box_upper < yop.box_constraint
    % Box constraint - upper bound
    
    properties (Constant)
        name = 'box_upper'
    end
    
    methods
        function obj = box_upper(var, bnd)
            obj@yop.box_constraint(var, bnd);
        end
    end
end