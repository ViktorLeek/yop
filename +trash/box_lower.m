classdef box_lower < yop.box_constraint
    % Box constraint - lower bound
    
    properties (Constant)
        name = 'box_lower'
    end
    
    methods
        function obj = box_lower(var, bnd)
            obj@yop.box_constraint(var, bnd);
        end
    end
end