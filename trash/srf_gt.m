classdef srf_gt < yop.srf
    
    properties (Constant)
        name = 'srf_gt'
    end
    
    methods
        function obj = srf_gt(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end
        
        function [type, constraint] = get_constraint(obj)
            if isa_variable(obj.lhs) && isnumeric(obj.rhs)
                % var >= num
                type = 'box';
                constraint = yop.box_lower(obj.lhs, obj.rhs);
                
            elseif isnumeric(obj.lhs) && isa_variable(obj.rhs)
                % num >= var
                % var <= num
                type = 'box';
                constraint = yop.box_upper(obj.rhs, obj.lhs);
                
            else
                type = 'inequality';
                constraint = [];
            end
        end
    end
end