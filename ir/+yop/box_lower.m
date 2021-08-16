classdef box_lower < yop.box_constraint
    % Box constraint - lower bound
    
    properties (Constant)
        name = 'box_lower'
    end
    
    methods
        function obj = box_lower(variable, bound)
            obj@yop.box_constraint(variable, bound);
        end
    end
end