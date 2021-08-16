classdef srf_ge < yop.srf
    
    properties (Constant)
        name = 'srf_ge'
    end
    
    methods
        function obj = srf_ge(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
    end
end