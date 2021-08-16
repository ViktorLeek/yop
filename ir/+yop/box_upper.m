classdef box_upper < yop.box_constraint
    % Box constraint - upper bound
    
    properties (Constant)
        name = 'box_upper'
    end
    
    methods
        function obj = box_upper(variable, bound)
            obj@yop.box_constraint(variable, bound);
        end
    end
end