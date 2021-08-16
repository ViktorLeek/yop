classdef srf_eq < yop.srf
    
    properties (Constant)
        name = 'srf_eq'
    end
    
    methods
        function obj = srf_eq(lhs, rhs)
            obj@yop.srf(lhs, rhs);
        end 
    end
end