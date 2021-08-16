classdef srf_gt < yop.srf
    
    properties (Constant)
        name = 'srf_gt'
    end
    
    methods
        function obj = srf_gt(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
    end
end