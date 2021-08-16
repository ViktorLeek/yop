classdef srf_le < yop.srf
    
    properties (Constant)
        name = 'srf_le'
    end
    
    methods
        function obj = srf_le(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
    end
end